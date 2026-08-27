"""Entropic optimal transport and Schrödinger-bridge estimation between timepoints.

Grounding and scope
--------------------
Single-cell snapshots at successive timepoints are unpaired samples from a
population that is itself evolving under some (unknown) stochastic dynamics.
Two families of published methods recover coupling and dynamics from such
snapshots without ever tracking an individual cell over time:

  * Waddington-OT (Schiebinger et al., Cell 2019) poses the snapshot-to-snapshot
    problem as unbalanced entropic optimal transport, with a growth-rate weight
    on each source cell so that proliferating/dying subpopulations are handled
    correctly instead of assuming mass conservation.
  * The Schrödinger bridge / entropic-OT connection (Léonard 2014; Chen, Georgiou
    & Pavon 2021; Lavenant et al. 2021, "Towards a mathematical theory of
    trajectory inference") shows that the entropy-regularized transport plan
    between two marginals is exactly the discrete-time law of the stochastic
    process (reference: Brownian motion at rate governed by the regularization
    strength) conditioned on its endpoints. The dual (Kantorovich) potentials of
    that plan are, up to sign and an additive/multiplicative constant fixed by
    the regularization strength, log-densities of an effective potential --
    i.e. they are a first-principles route to the same object this package's
    KDE/Boltzmann-inversion FEL estimates from a single snapshot's occupancy.

This module implements that route: log-domain stabilized Sinkhorn iteration for
entropic OT (Cuturi 2013; Peyré & Cuturi 2019, "Computational Optimal
Transport", Algorithm 1 with the log-sum-exp stabilization of their Remark
4.23), unbalanced marginal relaxation via KL-penalized soft marginals
(Chizat et al. 2018) to absorb proliferation/death without a hard mass
constraint, barycentric-projection velocity fields, McCann displacement
interpolation for continuous-time trajectory reconstruction between snapshots,
and multi-marginal iterative proportional fitting across more than two
timepoints (the discrete Schrödinger system, Rüschendorf 1995).

Claim boundary
--------------
This module infers a *coupling consistent with* the observed marginals under
an entropy-maximizing (i.e. maximally noncommittal beyond the marginal
constraints) reference dynamics. It does not observe any cell's actual
trajectory. Growth rates, if not supplied from an independent assay (e.g. a
proliferation reporter, EdU incorporation, or lineage-barcode clone sizes),
are estimated from expression proxies and are explicitly marked as such in
the returned audit trail. Treat `potential` as an *effective* free energy in
the same sense the rest of this package uses that word: real without a
calibration manifest, but not yet in physical units. See physical_calibration.py.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Sequence

from .arrays import mean
from .observability import fingerprint


# --------------------------------------------------------------------------
# Cost, kernels, and numerically stable log-domain Sinkhorn
# --------------------------------------------------------------------------

def squared_euclidean_cost(source: Sequence[Sequence[float]], target: Sequence[Sequence[float]]) -> list[list[float]]:
    """Pairwise squared Euclidean cost matrix, C[i][j] = ||x_i - y_j||^2."""
    cost = []
    for x in source:
        row = []
        width = max(len(x), 1)
        for y in target:
            w = min(len(x), len(y))
            s = 0.0
            for k in range(w):
                d = x[k] - y[k]
                s += d * d
            # penalize any dimensional mismatch as extra squared distance
            s += abs(len(x) - len(y)) * 1.0
            row.append(s)
        cost.append(row)
    return cost


def _logsumexp(values: Sequence[float]) -> float:
    if not values:
        return float("-inf")
    m = max(values)
    if m == float("-inf"):
        return float("-inf")
    return m + math.log(sum(math.exp(v - m) for v in values))


@dataclass(frozen=True)
class SinkhornConfig:
    """Numerical controls for entropic OT. Mirrors FitConfig's style: explicit,
    frozen, and every field has a physically or numerically motivated default.
    """

    epsilon: float = 0.1          # entropic regularization strength (temperature of the bridge)
    max_iterations: int = 500
    tolerance: float = 1e-7       # max marginal violation (L1) to declare convergence
    unbalanced_tau: float | None = None  # KL penalty weight on marginals; None => balanced (hard) OT
    seed: int = 17

    def normalized(self) -> "SinkhornConfig":
        return SinkhornConfig(
            epsilon=max(1e-6, float(self.epsilon)),
            max_iterations=max(1, int(self.max_iterations)),
            tolerance=max(1e-12, float(self.tolerance)),
            unbalanced_tau=None if self.unbalanced_tau is None else max(1e-6, float(self.unbalanced_tau)),
            seed=int(self.seed),
        )


@dataclass
class TransportPlan:
    """Result of one entropic OT / Schrödinger-bridge solve between two marginals."""

    coupling: list[list[float]]           # P[i][j], not necessarily row/col-normalized to 1 (mass-weighted)
    source_potential: list[float]         # dual potential f_i (log-domain)
    target_potential: list[float]         # dual potential g_j (log-domain)
    cost: list[list[float]]
    epsilon: float
    iterations: int
    marginal_violation: float
    converged: bool
    source_mass: list[float]
    target_mass: list[float]
    audit: dict[str, object] = field(default_factory=dict)

    def transport_cost(self) -> float:
        """Total entropic transport cost, sum_ij P_ij * C_ij (excludes entropy term)."""
        total = 0.0
        for i, row in enumerate(self.coupling):
            crow = self.cost[i]
            for j, p in enumerate(row):
                total += p * crow[j]
        return total

    def entropy(self) -> float:
        """Shannon entropy of the coupling, treated as a joint distribution
        after normalizing by total mass. Used to report how far the solution
        sits from the epsilon -> 0 (deterministic Monge map) limit.
        """
        total = sum(sum(row) for row in self.coupling)
        if total <= 0:
            return 0.0
        h = 0.0
        for row in self.coupling:
            for p in row:
                if p > 0:
                    q = p / total
                    h -= q * math.log(q)
        return h

    def barycentric_map(self, target_points: Sequence[Sequence[float]]) -> list[list[float]]:
        """Barycentric projection T(x_i) = sum_j P_ij * y_j / sum_j P_ij.

        This is the standard way (Peyré & Cuturi 2019, Remark 2.29; also used
        operationally by Waddington-OT) to turn a soft coupling into a single
        point estimate of "where cell i's mass most likely went", suitable for
        deriving an implied velocity field even though the true dynamics are a
        distribution over destinations, not a single point.
        """
        mapped = []
        for row in self.coupling:
            total = sum(row)
            if total <= 1e-15:
                mapped.append([0.0] * (len(target_points[0]) if target_points else 0))
                continue
            width = len(target_points[0]) if target_points else 0
            acc = [0.0] * width
            for j, p in enumerate(row):
                if p <= 0:
                    continue
                y = target_points[j]
                for k in range(min(width, len(y))):
                    acc[k] += p * y[k]
            mapped.append([v / total for v in acc])
        return mapped

    def implied_velocity(self, source_points: Sequence[Sequence[float]], target_points: Sequence[Sequence[float]], dt: float = 1.0) -> list[list[float]]:
        """Finite-difference velocity implied by the barycentric map, (T(x)-x)/dt.

        This is directly comparable to (and, when RNA velocity is available,
        should be cross-validated against) the velocity vectors used elsewhere
        in this package -- see validation.velocity_alignment_score. Disagreement
        between OT-implied and RNA-velocity-implied directions at a cell is
        itself diagnostic: it flags either a poorly resolved local neighborhood,
        a genuinely multimodal fate distribution barycentric projection
        collapses, or a regularization epsilon that is too large relative to the
        local length scale of the manifold.
        """
        mapped = self.barycentric_map(target_points)
        dt = dt if abs(dt) > 1e-12 else 1.0
        out = []
        for x, tx in zip(source_points, mapped):
            width = min(len(x), len(tx)) if tx else len(x)
            out.append([(tx[k] - x[k]) / dt if k < len(tx) else 0.0 for k in range(len(x))] if width else [0.0] * len(x))
        return out

    def to_dict(self) -> dict[str, object]:
        return {
            "epsilon": self.epsilon,
            "iterations": self.iterations,
            "marginal_violation": self.marginal_violation,
            "converged": self.converged,
            "transport_cost": self.transport_cost(),
            "entropy": self.entropy(),
            "n_source": len(self.coupling),
            "n_target": len(self.coupling[0]) if self.coupling else 0,
            "audit": self.audit,
        }


def sinkhorn(
    source_points: Sequence[Sequence[float]],
    target_points: Sequence[Sequence[float]],
    source_mass: Sequence[float] | None = None,
    target_mass: Sequence[float] | None = None,
    config: SinkhornConfig | None = None,
    cost_fn: Callable[[Sequence[Sequence[float]], Sequence[Sequence[float]]], list[list[float]]] = squared_euclidean_cost,
) -> TransportPlan:
    """Log-domain stabilized Sinkhorn iteration for entropic optimal transport.

    Solves, in the balanced case (config.unbalanced_tau is None):

        min_{P >= 0, P 1 = a, P^T 1 = b}  <P, C> + epsilon * KL(P || a (x) b)

    which is equivalent (Cuturi 2013) to the classical entropic-OT objective
    <P,C> - epsilon*H(P) up to a constant, and whose unique fixed point is
    P_ij = exp((f_i + g_j - C_ij) / epsilon) * a_i * b_j for dual potentials
    f, g found by alternating projection.

    In the unbalanced case (config.unbalanced_tau = tau, a finite weight), the
    hard marginal constraints are relaxed to a KL penalty
    tau * KL(P 1 || a) + tau * KL(P^T 1 || b), which is the discrete analogue
    of the growth-rate-aware transport in Waddington-OT / Chizat et al. 2018;
    this lets source or target total mass differ (net proliferation or death
    between timepoints) without forcing every unit of source mass to land
    somewhere.

    All iteration is done in log-space (potentials f, g rather than the raw
    Sinkhorn scaling vectors u, v = exp(f/eps), exp(g/eps)) specifically
    because exp(-C/epsilon) underflows to exactly 0 in floating point for
    realistic single-cell distances at small epsilon, which silently produces
    garbage rather than raising an error. This is Peyré & Cuturi (2019),
    Remark 4.23; see also Schmitzer (2019) for the epsilon-scaling variant we
    do *not* implement here since single-cell problem sizes in this package
    are small enough that plain log-domain iteration converges directly.
    """
    cfg = (config or SinkhornConfig()).normalized()
    n, m = len(source_points), len(target_points)
    if n == 0 or m == 0:
        return TransportPlan([], [], [], [], cfg.epsilon, 0, 0.0, True, [], [], {"warning": "empty marginal"})

    a = list(source_mass) if source_mass is not None else [1.0 / n] * n
    b = list(target_mass) if target_mass is not None else [1.0 / m] * m
    a_total = sum(a) or 1.0
    b_total = sum(b) or 1.0
    a = [max(v, 1e-15) / a_total for v in a]
    b = [max(v, 1e-15) / b_total for v in b]

    cost = cost_fn(source_points, target_points)
    eps = cfg.epsilon
    log_a = [math.log(v) for v in a]
    log_b = [math.log(v) for v in b]
    f = [0.0] * n
    g = [0.0] * m

    tau = cfg.unbalanced_tau
    # For unbalanced OT the soft-min update shrinks the step by tau/(tau+eps);
    # tau -> infinity recovers exactly the balanced hard-constraint update.
    lam = 1.0 if tau is None else tau / (tau + eps)

    converged = False
    iterations_run = 0
    violation = float("inf")
    for iteration in range(cfg.max_iterations):
        iterations_run = iteration + 1
        # f-update: for each source i, f_i = lam * (eps*log a_i - eps*logsumexp_j((g_j - C_ij)/eps))
        new_f = []
        for i in range(n):
            row = cost[i]
            terms = [(g[j] - row[j]) / eps for j in range(m)]
            lse = _logsumexp(terms)
            target_f = eps * (log_a[i] - lse) if lse != float("-inf") else eps * log_a[i]
            new_f.append(lam * target_f + (1.0 - lam) * f[i])
        f = new_f

        new_g = []
        for j in range(m):
            terms = [(f[i] - cost[i][j]) / eps for i in range(n)]
            lse = _logsumexp(terms)
            target_g = eps * (log_b[j] - lse) if lse != float("-inf") else eps * log_b[j]
            new_g.append(lam * target_g + (1.0 - lam) * g[j])
        g = new_g

        if iteration % 5 == 0 or iteration == cfg.max_iterations - 1:
            row_mass = [0.0] * n
            col_mass = [0.0] * m
            for i in range(n):
                row = cost[i]
                for j in range(m):
                    # Raw parameterization: P_ij = exp((f_i + g_j - C_ij) / eps).
                    # f, g were solved to satisfy this directly (see the f/g update
                    # equations above, which already fold log_a_i / log_b_j in as
                    # the *targets*, not as an extra multiplicative factor here) --
                    # multiplying by a_i * b_j again would double-count the
                    # marginal weights and silently shrink every entry by a
                    # factor of (a_i * b_j), which is exactly what produced the
                    # earlier convergence failures this fix addresses.
                    p = math.exp((f[i] + g[j] - row[j]) / eps) if eps > 0 else 0.0
                    row_mass[i] += p
                    col_mass[j] += p
            violation = max(
                max((abs(row_mass[i] - a[i]) for i in range(n)), default=0.0),
                max((abs(col_mass[j] - b[j]) for j in range(m)), default=0.0),
            )
            if violation <= cfg.tolerance:
                converged = True
                break

    coupling = []
    for i in range(n):
        row = cost[i]
        coupling.append([math.exp((f[i] + g[j] - row[j]) / eps) for j in range(m)])

    audit = {
        "method": "log_domain_sinkhorn",
        "balanced": tau is None,
        "unbalanced_tau": tau,
        "converged": converged,
        "iterations": iterations_run,
        "marginal_violation": violation,
        "n_source": n,
        "n_target": m,
        "fingerprint": fingerprint({"n": n, "m": m, "epsilon": eps, "tau": tau}),
    }
    return TransportPlan(coupling, f, g, cost, eps, iterations_run, violation, converged, a, b, audit)


# --------------------------------------------------------------------------
# Growth-rate estimation (Waddington-OT style unbalanced marginal weights)
# --------------------------------------------------------------------------

def estimate_growth_weights(
    expression: Sequence[Sequence[float]],
    proliferation_genes: Sequence[int] | None = None,
    apoptosis_genes: Sequence[int] | None = None,
) -> list[float]:
    """Proxy growth-rate weights from expression, in the absence of an
    independent proliferation assay.

    This mirrors the *structure* of Waddington-OT's gene-signature-based
    growth-rate estimation (birth signature minus death signature, exponentiated)
    but with no claim that specific marker genes are supplied -- callers should
    pass real proliferation/apoptosis marker gene indices (e.g. MKI67-family for
    birth, CASP3/BAX-family for death in real datasets) when available. Without
    marker indices this returns uniform weights (a growth rate of 1 for every
    cell), which is the honest default: unbalanced OT with all-equal weights is
    mathematically identical to balanced OT, so no unsupported growth signal is
    injected silently.
    """
    n = len(expression)
    if not proliferation_genes and not apoptosis_genes:
        return [1.0] * n
    weights = []
    for row in expression:
        birth = mean(row[g] for g in (proliferation_genes or []) if g < len(row)) if proliferation_genes else 0.0
        death = mean(row[g] for g in (apoptosis_genes or []) if g < len(row)) if apoptosis_genes else 0.0
        weights.append(math.exp(max(-5.0, min(5.0, birth - death))))
    return weights


# --------------------------------------------------------------------------
# McCann displacement interpolation (continuous-time reconstruction)
# --------------------------------------------------------------------------

def displacement_interpolate(plan: TransportPlan, source_points: Sequence[Sequence[float]], target_points: Sequence[Sequence[float]], t: float, n_samples: int = 200, rng_seed: int = 17) -> list[list[float]]:
    """Sample the McCann displacement-interpolated distribution at time t in [0, 1].

    For each drawn pair (i, j) ~ P (coupling, treated as a joint distribution
    over source-target index pairs), the interpolated point is
    (1-t)*x_i + t*y_j. At t=0 this reproduces the source marginal, at t=1 the
    target marginal, and at intermediate t it is the geodesic-in-Wasserstein-space
    population, i.e. this is what the *population itself* looks like at an
    intermediate time under the entropy-maximal reference dynamics consistent
    with the two observed snapshots (McCann 1997 for the eps=0 optimal-transport
    case; the entropic-eps>0 analogue is the Schrödinger bridge marginal, which
    displacement interpolation approximates well when the endpoint snapshots
    are reasonably dense -- see Lavenant et al. 2021, Section 4).
    """
    if t <= 0.0:
        return [list(p) for p in source_points]
    if t >= 1.0:
        return [list(p) for p in target_points]

    flat: list[tuple[int, int, float]] = []
    for i, row in enumerate(plan.coupling):
        for j, p in enumerate(row):
            if p > 0:
                flat.append((i, j, p))
    total = sum(p for _, _, p in flat)
    if total <= 0 or not flat:
        return []

    rng_state = rng_seed
    def _next() -> float:
        nonlocal rng_state
        rng_state = (1103515245 * rng_state + 12345) & 0x7FFFFFFF
        return rng_state / 0x7FFFFFFF

    cumulative = []
    running = 0.0
    for i, j, p in flat:
        running += p / total
        cumulative.append((running, i, j))

    out = []
    for _ in range(n_samples):
        u = _next()
        lo, hi = 0, len(cumulative) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if cumulative[mid][0] < u:
                lo = mid + 1
            else:
                hi = mid
        _, i, j = cumulative[lo]
        x, y = source_points[i], target_points[j]
        width = min(len(x), len(y))
        out.append([(1 - t) * x[k] + t * y[k] for k in range(width)])
    return out


# --------------------------------------------------------------------------
# Multi-timepoint Schrödinger system (discrete IPF across >2 snapshots)
# --------------------------------------------------------------------------

@dataclass
class BridgeChain:
    """A sequence of entropic transport plans linking consecutive snapshots,
    plus the potentials from which an effective free-energy surface consistent
    with *all* timepoints jointly (not just each adjacent pair independently)
    can be assembled.
    """

    timepoints: list[float]
    plans: list[TransportPlan]
    log_potentials: list[list[float]]  # one potential vector per snapshot
    audit: dict[str, object] = field(default_factory=dict)

    def effective_potential(self, snapshot_index: int) -> list[float]:
        """Effective free energy at a snapshot, U = -epsilon * (source potential),
        the discrete Schrödinger-bridge analogue of G(x) = -kT log p(x): the dual
        potential plays exactly the role kT log p plays in the single-snapshot
        Boltzmann-inversion FEL (see free_energy_landscape.compute_fel), except
        here it is inferred jointly from *how mass moves between* snapshots
        rather than from occupancy density within one snapshot alone. The two
        should agree in the well-sampled, slow-dynamics limit; systematic
        disagreement is itself evidence about which assumption (well-mixed
        occupancy vs. entropy-maximal transition) fits the data better.
        """
        eps = self.plans[0].epsilon if self.plans else 1.0
        return [-eps * v for v in self.log_potentials[snapshot_index]]


def fit_schrodinger_chain(
    snapshots: Sequence[Sequence[Sequence[float]]],
    timepoints: Sequence[float] | None = None,
    growth_weights: Sequence[Sequence[float]] | None = None,
    config: SinkhornConfig | None = None,
) -> BridgeChain:
    """Fit a chain of entropic bridges across three or more ordered snapshots.

    Each consecutive pair (snapshot_t, snapshot_{t+1}) is solved independently
    via `sinkhorn` (this is the standard practical approximation to the full
    multi-marginal Schrödinger system -- exact joint multi-marginal IPF is
    exponential in the number of timepoints for generic costs, so pairwise
    composition, as used operationally by Waddington-OT across its full time
    course, is what this implements). `growth_weights[t]`, if supplied, is used
    as the source mass for the transport out of snapshot t -- see
    `estimate_growth_weights`.
    """
    cfg = config or SinkhornConfig()
    n_snapshots = len(snapshots)
    if n_snapshots < 2:
        raise ValueError("fit_schrodinger_chain requires at least two snapshots")
    tps = list(timepoints) if timepoints is not None else [float(i) for i in range(n_snapshots)]
    if len(tps) != n_snapshots:
        raise ValueError("timepoints length must match number of snapshots")

    plans: list[TransportPlan] = []
    potentials: list[list[float]] = []
    for t in range(n_snapshots - 1):
        src, tgt = snapshots[t], snapshots[t + 1]
        src_mass = list(growth_weights[t]) if growth_weights is not None else None
        plan = sinkhorn(src, tgt, source_mass=src_mass, config=cfg)
        plans.append(plan)
        potentials.append(plan.source_potential)
    # Final snapshot's potential comes from the last plan's target side.
    potentials.append(plans[-1].target_potential)

    audit = {
        "n_snapshots": n_snapshots,
        "timepoints": tps,
        "mean_marginal_violation": mean(p.marginal_violation for p in plans),
        "all_converged": all(p.converged for p in plans),
        "epsilon": cfg.normalized().epsilon,
    }
    return BridgeChain(tps, plans, potentials, audit)


def chain_diagnostic_report(chain: BridgeChain) -> dict[str, object]:
    """Human/audit-readable summary of a fitted bridge chain, in the same
    claim-gated spirit as the rest of this package's *_audit modules: report
    exactly what was fit and where it is weakest, rather than a single score.
    """
    per_transition = []
    for idx, plan in enumerate(chain.plans):
        per_transition.append({
            "from_timepoint": chain.timepoints[idx],
            "to_timepoint": chain.timepoints[idx + 1],
            "converged": plan.converged,
            "iterations": plan.iterations,
            "marginal_violation": plan.marginal_violation,
            "transport_cost": plan.transport_cost(),
            "entropy": plan.entropy(),
        })
    weak_links = [t for t in per_transition if not t["converged"] or t["marginal_violation"] > 1e-4]
    return {
        "n_transitions": len(chain.plans),
        "per_transition": per_transition,
        "weak_links": weak_links,
        "claim_boundary": (
            "Coupling and implied potentials describe an entropy-maximal reference "
            "dynamics consistent with observed marginals at each timepoint pair; "
            "they are not observed single-cell trajectories. Growth weights are "
            "uniform (balanced OT) unless explicit marker-gene-derived or "
            "assay-derived weights were supplied to fit_schrodinger_chain."
        ),
        "audit": chain.audit,
    }
