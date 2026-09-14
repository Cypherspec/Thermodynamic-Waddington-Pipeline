# Preregistration: a falsifiable test of the commitment prediction

This document states, in advance, a specific prediction from the Thermodynamic
Waddington pipeline and how an experiment would confirm or falsify it. It is
written before any perturbation data is collected. It is deliberately built on
the pipeline's robust, cross-dataset result (entropy production and the
commitment coordinate), not on the dataset-dependent ordering result.

## What the pipeline claims, and how well supported it is

- **Supported across two datasets:** the differentiation process is
  measurably irreversible. Entropy production separates directed dynamics from a
  velocity-shuffled control on pancreatic endocrinogenesis (directed p ~ 0.01)
  and on gastrulation erythroid (directed p ~ 0.01), while the shuffled control
  is not significant on either.
- **Supported on pancreas:** the committor places fate commitment at the
  Pre-endocrine stage (committor crosses 0.5), with a potential-of-mean-force
  barrier of about 1.6 kT [1.4, 1.9] along the reaction coordinate.
- **Not claimed:** that the committor orders cells better than PCA-based
  pseudotime in general. It does on pancreas and does not on gastrulation
  erythroid. The prediction below does not depend on ordering.

## Primary hypothesis (H1)

The commitment barrier sits at the transition state (committor q ~ 0.5,
Pre-endocrine on pancreas). Therefore a perturbation applied to cells **at the
transition state** will change the terminal-fate outcome **more** than the same
perturbation applied to cells **at the progenitor state** (q ~ 0, Ductal).

This is the operational meaning of a barrier: the system is most sensitive where
the free energy is highest along the path.

## Pre-specified design

- **System:** an in vitro pancreatic endocrine differentiation with a Beta-fate
  readout, or an equivalent lineage with a defined progenitor, transition, and
  terminal population.
- **Arms:** (A) perturb at progenitor state, (B) perturb at transition state,
  (C) vehicle control. The perturbation is a defined, titratable push on the
  fate axis (a small-molecule or genetic modulator chosen independently of this
  pipeline).
- **Readout:** terminal-fate fraction (e.g. Beta+ cells) by scRNA-seq or flow,
  measured at a fixed endpoint.
- **Replication:** at least 3 independent donors or biological replicates per
  arm, sample sizes fixed in advance by a power analysis for the effect size
  below.

## Pre-specified analysis and decision rule

- **Effect measure:** change in terminal-fate fraction versus control, per arm.
- **Confirmation of H1:** arm B (transition-state perturbation) shows a larger
  absolute change than arm A (progenitor perturbation), with a pre-registered
  minimum difference of 10 percentage points and a paired test at alpha = 0.05,
  donor-paired.
- **Falsification of H1:** arm A greater than or equal to arm B, or neither arm
  differs from control. Either outcome falsifies the barrier-location claim and
  is to be reported as such.
- No optional stopping, no post-hoc arm selection, no swapping the readout.

## What a null result would mean

A null or reversed result would show that the pipeline's committor barrier does
not locate the point of maximal fate sensitivity. That is a genuine negative and
would be reported, in the same spirit as the documented dentate gyrus
non-replication and the gastrulation ordering result.

## Analysis code

The pipeline functions used to generate the prediction, fixed at the current
commit: `developmental_coordinate`, `commitment_profile`,
`committor_free_energy_profile`, and `estimate_entropy_production`. The
prediction (Pre-endocrine transition state, ~1.6 kT barrier) is reproducible via
`python benchmarks/commitment_barrier.py`.
