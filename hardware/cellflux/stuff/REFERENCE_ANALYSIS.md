# Reference analysis implementation — NOT EXECUTED

This is source text for a proposed small-state continuous-time Markov analysis module, not a deployed Python package or the supplied Thermodynamic Waddington engine. No Python execution tool is available in this CAD workspace. The assertions below are analytic test fixtures, not passed tests. Extract and independently review/test before use. Dense linear algebra is intentional for small state graphs; this is not a single-cell-scale production solver.

## Input contract
`L[i,j]` is the transition rate from i to j, rows sum to zero. Its clock must come from calibrated longitudinal observations or a separately validated kinetic model. A row-stochastic RNA-velocity transition matrix is NOT this input and must not be silently interpreted as rates per second. A and B are nonempty, disjoint, prespecified state-index sets. This reference deliberately rejects reducible models: absorbing developmental systems require a different nonstationary analysis, not artificial reverse edges.

```python
import numpy as np


def checked_solve(a, b):
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("nonfinite linear system")
    condition = np.linalg.cond(a)
    if not np.isfinite(condition) or condition > 1e12:
        raise ValueError("ill-conditioned system; do not publish estimates")
    x = np.linalg.solve(a, b)
    if not np.isfinite(x).all() or not np.allclose(a @ x, b, rtol=1e-9, atol=1e-10):
        raise ValueError("linear solve residual failed")
    return x


def checked_generator(rates):
    l = np.array(rates, dtype=float, copy=True)
    if l.ndim != 2 or l.shape[0] != l.shape[1] or l.shape[0] < 2:
        raise ValueError("generator must be square with at least two states")
    if not np.isfinite(l).all():
        raise ValueError("nonfinite rates")
    off = l.copy()
    np.fill_diagonal(off, 0.0)
    if (off < 0).any() or (np.diag(l) > 0).any():
        raise ValueError("invalid rate signs")
    scale = float(np.max(-np.diag(l)))
    if scale <= 0 or not np.allclose(l.sum(axis=1) / scale, 0, atol=1e-12, rtol=0):
        raise ValueError("invalid row sums or zero generator")
    adjacency = off > 0
    for graph in (adjacency, adjacency.T):
        seen, stack = {0}, [0]
        while stack:
            i = stack.pop()
            for j in np.flatnonzero(graph[i]):
                j = int(j)
                if j not in seen:
                    seen.add(j)
                    stack.append(j)
        if len(seen) != l.shape[0]:
            raise ValueError("reducible model: stationary reference is inapplicable")
    return l, scale


def boundary_indices(values, n):
    values = list(values)
    if not values or any(isinstance(v, (bool, np.bool_)) or not isinstance(v, (int, np.integer)) for v in values):
        raise ValueError("boundary must contain integer state indices")
    out = np.unique(values)
    if out.min() < 0 or out.max() >= n:
        raise ValueError("boundary index outside state space")
    return out


def dirichlet(generator, zero, one):
    n = len(generator)
    interior = np.setdiff1d(np.arange(n), np.union1d(zero, one))
    q = np.zeros(n)
    q[one] = 1.0
    if len(interior):
        q[interior] = checked_solve(
            generator[np.ix_(interior, interior)],
            -generator[np.ix_(interior, one)].sum(axis=1),
        )
    if (q < -1e-9).any() or (q > 1 + 1e-9).any():
        raise ValueError("committor out of bounds")
    return np.clip(q, 0.0, 1.0)


def analyze_ctmc(rates, a_states, b_states, *, time_unit="model_time"):
    if not isinstance(time_unit, str) or not time_unit.strip():
        raise ValueError("explicit time unit required")
    l, scale = checked_generator(rates)
    n = len(l)
    a = boundary_indices(a_states, n)
    b = boundary_indices(b_states, n)
    if np.intersect1d(a, b).size:
        raise ValueError("A and B overlap")
    g = l / scale
    system = g.T.copy()
    system[-1, :] = 1.0
    rhs = np.zeros(n)
    rhs[-1] = 1.0
    pi = checked_solve(system, rhs)
    if (pi <= 0).any() or not np.allclose(pi @ g, 0, atol=1e-10, rtol=0):
        raise ValueError("invalid stationary distribution")
    pi /= pi.sum()
    q_plus = dirichlet(g, a, b)
    # Time reversal: Lrev_ij = pi_j L_ji / pi_i.
    reverse = (g.T * pi[None, :]) / pi[:, None]
    q_minus = dirichlet(reverse, b, a)
    outside_b = np.setdiff1d(np.arange(n), b)
    mfpt = np.zeros(n)
    mfpt[outside_b] = checked_solve(
        g[np.ix_(outside_b, outside_b)], -np.ones(len(outside_b))
    ) / scale
    if not np.isfinite(mfpt).all() or (mfpt < -1e-10).any():
        raise ValueError("invalid first-passage times")
    directed_flux = pi[:, None] * l
    np.fill_diagonal(directed_flux, 0.0)
    current = directed_flux - directed_flux.T
    entropy_rate = 0.0
    missing_reverse_edges = []
    for i in range(n):
        for j in range(i + 1, n):
            forward, backward = directed_flux[i, j], directed_flux[j, i]
            if forward == 0 and backward == 0:
                continue
            if forward == 0 or backward == 0:
                missing_reverse_edges.append((i, j))
                entropy_rate = float("inf")
            elif np.isfinite(entropy_rate):
                entropy_rate += (forward - backward) * (np.log(forward) - np.log(backward))
    reactive = q_minus[:, None] * directed_flux * q_plus[None, :]
    reactive_net = np.maximum(reactive - reactive.T, 0.0)
    interior = np.setdiff1d(np.arange(n), np.union1d(a, b))
    if not np.allclose((reactive.sum(axis=1) - reactive.sum(axis=0))[interior] / scale, 0, atol=1e-9):
        raise ValueError("reactive flux conservation failed")
    dimensionless_potential = -np.log(pi)
    dimensionless_potential -= dimensionless_potential.min()
    return {
        "stationary_probability": pi,
        "dimensionless_potential": dimensionless_potential,
        "forward_committor": q_plus,
        "backward_committor": q_minus,
        "mfpt_to_b": mfpt,
        "stationary_directed_flux": directed_flux,
        "stationary_current": current,
        "coarse_grained_entropy_rate": entropy_rate,
        "missing_reverse_edges": missing_reverse_edges,
        "reactive_gross_flux": reactive,
        "reactive_net_flux": reactive_net,
        "reactive_events_per_time": float(reactive[a, :].sum()),
        "time_unit": time_unit,
        "entropy_rate_unit": "nats/" + time_unit,
        "validation_status": "REFERENCE_IMPLEMENTATION_NOT_BIOLOGICALLY_VALIDATED",
    }


def reference_fixtures():
    # These assertions have NOT been run in this workspace.
    reversible = [[-1, 1, 0], [1, -2, 1], [0, 1, -1]]
    r = analyze_ctmc(reversible, [0], [2])
    np.testing.assert_allclose(r["stationary_probability"], [1/3]*3)
    np.testing.assert_allclose(r["forward_committor"], [0, 0.5, 1])
    np.testing.assert_allclose(r["backward_committor"], [1, 0.5, 0])
    np.testing.assert_allclose(r["mfpt_to_b"], [3, 2, 0])
    np.testing.assert_allclose(r["coarse_grained_entropy_rate"], 0, atol=1e-12)
    np.testing.assert_allclose(r["reactive_events_per_time"], 1/6)
    driven = [[-3, 2, 1], [1, -3, 2], [2, 1, -3]]
    r = analyze_ctmc(driven, [0], [2])
    np.testing.assert_allclose(r["stationary_probability"], [1/3]*3)
    np.testing.assert_allclose(r["forward_committor"], [0, 2/3, 1])
    np.testing.assert_allclose(r["backward_committor"], [1, 2/3, 0])
    np.testing.assert_allclose(r["mfpt_to_b"], [5/7, 4/7, 0])
    np.testing.assert_allclose(r["coarse_grained_entropy_rate"], np.log(2))
    for invalid in ([[0, 0], [0, 0]], [[0, 0], [1, -1]]):
        try:
            analyze_ctmc(invalid, [0], [1])
        except ValueError:
            pass
        else:
            raise AssertionError("invalid/reducible generator accepted")
```

## Interpretation and deployment gates
The potential is dimensionless negative log stationary state mass, depends on state discretization, and is not physical cellular free energy. Entropy rate uses positive directed fluxes, never the logarithm of signed net current. A missing observed reverse transition yields an infinite estimator with a diagnostic, not evidence of infinite cellular dissipation. Bootstrap whole independent lineages/experiments and vary state definition, lag and kinetic estimation; do not bootstrap correlated frames as independent cells. Inferred values apply to this model and its assumptions, not direct calorimetry.

Reactive event frequency is not automatically inverse MFPT. Nonreversible backward committors require the time-reversed process. Do not replace them with one minus forward committor. Cycle decomposition, rate estimation, RNA preprocessing, tracking, uncertainty quantification, hardware drivers, UI and original-engine adapter remain separate unimplemented modules.

Grounding: deeptime authors' Transition path theory documentation (boundary-value and reactive-flux definitions); NumPy `numpy.linalg.solve` documentation; Ghosal and Bisker, arXiv:2205.14688 (partial observation/coarse-graining limits). These sources do not validate this unexecuted code.
