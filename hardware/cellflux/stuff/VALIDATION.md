# Verification and release gates — Revision B

**Status: executable concept CAD; not a working, certified or fabrication-released instrument.** Proposed tests below have NOT been performed. Acceptance numbers are candidate requirements to approve against the chosen cell type and assay, not measured performance.

| Gate | Evidence / proposed acceptance | Current status |
|---|---|---|
| Source integrity | Preserve original assembly bodies; separate files/parts; no final assembly fusion | Retained; revision B adds interfaces/routes |
| Solver constraints | No free or conflicting sketch variables in edited/new part entrypoints | Checked by CAD tools; main imports report zero local sketches |
| Geometric execution | Main and new/edited part entrypoints execute | CAD execution checked; not proof of fit/function |
| Visual features | Selected channel/fiducial voids, open tube lumen, lid penetrations, service routing, restored covers | Targeted views reviewed; not exhaustive collision analysis |
| Rate-model mathematics | Execute REFERENCE_ANALYSIS.md fixtures plus independent solver comparison and randomized/property tests | Reference text only, unexecuted |
| Original-engine integration | Freeze actual API, units, model versions, input fixtures and expected outputs | Engine source/API not supplied |
| Identity and provenance | No silent identity reassignment; reject ambiguous endpoint registration; immutable frame/event/calibration links | IMPLEMENTATION.md specification only |
| Time synchronization | Candidate exposure-to-event clock error <=1 ms, measured over a full run; characterize fluid arrival separately | No hardware or measured timing |
| Registration | Candidate held-out registration error <=2 um across field after cartridge reseating; do not assess only the three fitted fiducials | No optical tests |
| Cell tracking | Manually annotated held-out lineage data; report identity switches, missed divisions and censoring; threshold agreed before test | No dataset or tracker implemented |
| Environmental control | Candidate 37.0 +/-0.5 C at sample for a mammalian assay; CO2 setpoint chosen from medium requirements; map gradients and recovery | Thermal/gas components are envelopes, not validated systems |
| Flow | Calibrated liquid-specific sensor range; candidate <=5% steady flow error; bubble/clog tests and measured perturbation arrival distribution | No hydraulic simulation or flow bench test |
| Pressure protection | Rated component selection, relief and independent shutdown; approved proof/leak procedure | No pressure rating assigned; operation disabled by specification |
| Motion and optics | Vendor-defined XY/focus travel, working distance, collision envelopes, optical prescription and calibration | Incomplete motion hardware and optical envelopes |
| Cartridge fabrication | Dimensional/channel-depth metrology, bonded-interface validation, wetting and extractables review | Geometry only; fabrication untested |
| Biological assay | Approved cell-specific protocol, capture/doublet/viability criteria, endpoint linkage and controls | No cells tested; 48-chamber layout remains a candidate |
| Scientific validity | Held-out lineage prediction, calibrated time units, model adequacy, uncertainty and coarse-graining sensitivity | No experimental inference or validated performance |
| Safety | Appropriate external electrical/optical/pressure/containment and ethics review; records signed by responsible personnel | No compliance claim |

## Bench sequence, once engineering and approvals are complete
1. Unpowered mechanical fit inspection and independent interference/tolerance review.
2. Electrical/pressure/optical protection verification using approved methods before wet operation.
3. Compatible nonbiological tracer tests: leaks, wetting, bubbles, flow calibration, selector carryover and arrival latency.
4. Optical bead/registration/focus characterization and full-duration thermal/gas stability.
5. Hardware-in-loop injected faults: camera loss, stale sensors, disk-full, clock reset, valve failure, thermal excursion and power loss. Record detection latency, physical safe state and recovery authorization.
6. Approved cell pilot, lineage annotation and matched/endpoint molecular measurements. No claim of repeated destructive sequencing of the same live cell.
7. Independent scientific benchmark suite followed by blinded biological validation and documented release review.

## Release blockers, not waived by successful rendering
Pressure supply/relief plumbing, gas delivery/exhaust, sensor insertion, electrical harnesses/interlocks, fully defined XY/focus motion, vendor-specific optics, structural retention and full tolerance/interference analysis remain incomplete. Liquid hoses currently land on geometric interface envelopes; the selector's internal switching mechanism is not modeled. Direct hose-to-tail joins require selected connection hardware or a validated continuous-tube replacement. Grommets and window seals are nominal geometry, not proven seals. No operational settings or biologically meaningful thermodynamic outputs may be inferred from CAD validation alone.
