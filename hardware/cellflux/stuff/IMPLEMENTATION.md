# Instrument implementation specification — Revision B

Status: proposed engineering specification, not deployed software, a validated assay, or a safety release. The existing Thermodynamic Waddington source/API has not been supplied; integration cannot be tested against it. KCL executes CAD only. The available file tools write KCL and Markdown, not a runnable Python/C++ repository.

## Chosen assay architecture
1. Live imaging cartridge: continuous observation with immutable chamber and lineage IDs; endpoint molecular imaging in place is the proposed same-cell link. This requires experimentally validated fixation/retention and registration, not just matching locations.
2. Droplet cartridge: independent matched-population sequencing branch. Never label these observations as the same tracked cells without a separately validated physical barcode/recovery method.
3. Endpoint FISH is a targeted measurement, not automatically a whole-transcriptome or RNA-velocity dataset. Record panel, probes, detection limits and missing genes. Spliced/unspliced sequencing is an offline input to a separately validated inference workflow.

## Identity contract
- `experimentId`: immutable UUID, not a filename.
- `cartridgeId`: scan before priming; store cartridge geometry revision and assay lot.
- `chamberId`: C001 through C048 ordered in increasing chip-local X; nominal center X = -2.35 + 0.1*(index-1) mm, Y = 0.0895 mm.
- `trackId`, `parentTrackId`: independent of chamber occupancy. A chamber can contain multiple descendants; chamber ID alone is not cell identity.
- Every frame: exposure midpoint monotonic time, UTC time, camera frame sequence, pixel scale, stage encoder coordinates, illumination channel, exposure, objective, focal position, calibration hash.
- Every perturbation: command time, measured arrival time, valve identity, pressure/flow trace, reagent lot and concentration reference. Pump command time is not automatically cellular exposure time.
- Endpoint mapping: track ID, segmentation ID, registration transform, fiducial residual, ambiguity flag, operator review and molecular observation ID. Ambiguous/censored cells retain data but are excluded from same-cell claims.
- Sequencing records: physical sample ID, batch, collection interval, FASTQ provenance, count-layer names, normalization/model versions and relation=`matchedPopulation` unless proven otherwise.

## Registration interface
Chip-local underside fiducials are non-collinear and have different diameters: F1=(-3,-0.6), D=0.08 mm; F2=(3,-0.6), D=0.12 mm; F3=(-2.5,0.6), D=0.16 mm. Their CAD presence does not prove optical contrast. Use a calibrated handedness-aware chip-to-camera transform. Determine mirror orientation from the asymmetry. Three correspondences determine an affine transform but provide no independent validation; use additional cell-free landmarks or calibration targets to estimate residuals. Repeatedly image fiducials, do not infer stage repeatability from commanded positions.

## Proposed control state machine
`BOOT -> SELF_TEST -> IDLE -> PRIME -> LOAD -> TRACK -> ENDPOINT -> COMPLETE`
Any unsafe condition -> latched `FAULT`; no automatic restart. `TRACK -> PAUSE` is permitted only under a defined validated holding protocol. Returning from PAUSE starts a new acquisition segment and preserves the gap.
- BOOT: outputs disabled; verify configuration signature, device identities and storage availability.
- SELF_TEST: read actual switch, pressure, temperature and encoder values; no substitution of a successful API call for measured readiness.
- PRIME/LOAD: require approved recipe and containment state. Do not start with unspecified pressure/flow/temperature limits.
- TRACK: hardware-trigger acquisition; log measured flow and environmental conditions; detect dropped frames, tracking ambiguity and occlusions.
- ENDPOINT: explicit confirmation; irreversible assay transition. Stop claims of live observation, record time and lock prior images. Image registration precedes molecular assignment.
- COMPLETE: reconcile frames and samples, write checksums, audit manifest and immutable provenance.

## Safety and fault I/O contract
Proposed inputs: both door switches, emergency-stop feedback, pressure high limit, pressure transducers, waste-high switch, tray leak detector, chip and enclosure temperatures, gas-flow feedback, stage end switches, drive faults, light-source ready/fault, storage health.
Proposed outputs: pressure-enable, normally closed liquid isolation valves, controlled depressurization command, heater enable, source enable/shutter, camera trigger, stage commands and alarm.
Hardware requirements: emergency stop and hazardous-light/pressure isolation must not depend on this Python process, USB responsiveness or the GUI. Over-temperature protection must be independent of normal software control. Stale/disconnected safety telemetry is a fault. Depressurization must use a reviewed contained path, not uncontrolled venting of sample. Exact wiring, safe states and response times require component selection and hazard analysis.

## Services and interfaces
| Service | Input | Output / boundary |
|---|---|---|
| Device adapters | Device-specific SDKs and approved configuration | Typed measured values, acknowledgments, errors; no silent simulation fallback |
| Acquisition coordinator | Recipe, interlocks, hardware clocks | Trigger sequence, append-only events, image manifests |
| Tracking | Raw images and calibration | Segmentation, track/lineage graphs, uncertainty and censoring |
| Endpoint registration | Fiducials and live/fixed images | Validated correspondence table with rejects |
| Molecular import | FASTQ/counts or FISH spots | Separate assay-specific observations and QC |
| State model | Time-labelled observations and linkage quality | Generator/transition model, uncertainty, calibrated time units |
| Waddington adapter | Validated model plus A/B definitions | Metrics with units, assumptions, uncertainty and model hash |
| Operator GUI | Read-only measurements and state | Deliberate recipe requests; no direct actuator bypass |

Keep control and scientific inference separate: analysis failure must not disable safety. Use hardware timing for acquisition; network timestamps alone are insufficient for demonstrated synchronization. Capture actual dropped-frame counts and calibration residuals. Secrets, human sample identifiers and access controls belong in a reviewed deployment design.

## Scientific acceptance contract
- Define the state representation and observation process before inferring rates. Do not call -log density in an arbitrary UMAP embedding a measured physical free energy.
- For the reference finite-state model, report -log stationary probability as a dimensionless landscape. Rate units must come from calibrated elapsed time, not latent time alone.
- Stationary distributions need not imply detailed balance; a stationary driven cycle can have nonzero currents.
- Use positive directed stationary fluxes in logarithmic entropy-production expressions, not signed net currents. No silent pseudocounts to manufacture reverse transitions.
- Treat an absorbing differentiation dataset separately from an ergodic stationary model. The reference solver deliberately rejects reducible models rather than returning spurious stationary landscapes.
- Committors require explicit A and B definitions. MFPTs require a time-calibrated transition model and reachable targets. Report censoring, model uncertainty and held-out calibration.
- Bootstrap by independent experiment/lineage, not by treating adjacent frames or sibling cells as independent samples. Assess cell division and selective loss biases.
- No quantitative cell-fate percentage is released without prospective/held-out calibration. No thermodynamic entropy or power units are claimed from an unvalidated coarse-grained model.

## Sources informing the design boundary
- Lane et al., Cell Systems (2017), DOI 10.1016/j.cels.2017.03.010: integration of live signaling and same-cell RNA sequencing requires a specific measurement/linkage method.
- Bergen et al., Nature Biotechnology (2020), DOI 10.1038/s41587-020-0591-3: RNA velocity uses a transcription/splicing dynamics model.
- Metzner, Schuette and Vanden-Eijnden, SIAM Multiscale Modeling & Simulation (2009), DOI 10.1137/070699500: finite-state transition path theory.
- Ghosal and Bisker, arXiv:2205.14688: partial observations can provide lower bounds on total entropy production under the stated model assumptions.

These sources justify methodological distinctions, not the operation or performance of this CAD design.
