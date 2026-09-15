# Single-cell thermodynamics instrument — CAD handoff

## Status and scope
Revision B is an expanded **layout/concept CAD assembly**, not a fabrication release or a functioning instrument. It contains 88 separate part/module definitions, shared `geometry.kcl` and `fluidRouting.kcl` construction libraries, and `main.kcl` assembly placement. Commercial devices remain generic installation envelopes. No vendor interface verification, optical prescription, pressure rating, electrical design or certification is implied.

The supplied `deep-research-report.md` is the design brief, not validated engineering evidence. Its scientific, product-specification, biosafety and compliance claims are not adopted as verified facts.

## Layout and editable dimensions
- Base: 900 x 650 x 20 mm at Z=30; roof top Z=654 mm. Operator hardware projects outside the base footprint.
- Z up, front negative Y. Wet bay forward/right, dry electronics aft/right, microscopy left.
- Two 75 x 25 x 3 mm chip bodies on separate 0.17 mm coverglass, stepped carriers, clamps and eight retention screw/nut pairs.
- Imaging axis X=-230, Y=-65 mm. Second cartridge is 130 mm to the right; no validated motion covering both stations is claimed.
- Droplet candidate: 30 um channel depth/nozzle width, five 1.2 mm ports.
- Imaging candidate: 30 um depth, 48 dead-end chambers, 20 um chamber width, 100 um pitch, two ports. Capture performance untested.
- Three noncollinear underside imaging-chip fiducials: (-3,-0.6), (3,-0.6), (-2.5,0.6) mm; diameters 80,120,160 um; recess depth15 um. Intended for endpoint registration, not automatic proof of cell identity.
- Six reagent reservoirs, sample/waste vessels, eight valve and eight sensor installation envelopes.
- Nine 1 mm OD / 0.5 mm ID curved liquid routes have 10 mm nominal centerline bend radii. Tube ends connect geometrically to cartridge tails, selector port recesses and vessel inlets. Actual connection hardware, strain relief and bend specifications remain to be selected.
- Live feed crosses above condenser at Z=486.87 mm, using a separate3 mm clearance passage in its support arm. This avoids its original conflict with the illumination assembly.
- Chamber frame seals, optical window/annular seal/retainer, seven lid tube inserts, six reservoir and two collection inserts are nominal CAD geometry, not tested seals.
- Four PCB standoffs and eight backing nut envelopes added. Cosmetic thread intent only.
- `serviceLift`/`serviceSpread` displace covers for inspection only; restore both to zero for assembled configuration. They are not kinematic mates.
- Repeated dimension-driven profiles and hose routes use custom functions. Source dimensions are editable; sketches inside functions are not point-and-click editable.

## Liquid interface map
| Route | From | To |
|---|---|---|
| Live feed | Selector common (90,-65) | Live chip inlet (-260,-65) |
| Live drain | Live chip outlet (-200,-65) | Waste lid (197,-210) |
| Cell feed | Front-left reservoir liquid port (193,-140) | Droplet cell inlet (-130,-65) |
| Front oil | Front-middle reservoir liquid port (263,-140) | Oil inlet (-100,-73) |
| Rear oil | Rear-middle reservoir liquid port (263,-70) | Oil inlet (-100,-57) |
| Bead feed | Front-right reservoir liquid port (333,-140) | Bead inlet (-115,-57) |
| Droplet collection | Droplet outlet (-70,-65) | Sample lid (313,-210) |
| Media | Rear-right reservoir liquid port (333,-70) | Selector input (70,-85) |
| Perturbation | Rear-left reservoir liquid port (193,-70) | Selector input (70,-45) |

Coordinates in mm. Selector is a purchased-module envelope with isolated port recesses, NOT a functioning modeled valve. Other four selector inputs are plugged. Reservoir pressure ports and collection vents still require rated plumbing/filters. Live-chip effluent goes to waste, not a lineage-resolved cell sorter. Endpoint imaging/FISH registration and matched-population sequencing are the intended data-linkage options; same-cell sequencing linkage is not implemented.

## Validation boundary
Revision A documented successful checks of its67 part entrypoints. Revision B new/edited part entrypoints were checked for fully constrained sketches and execution; main's check reports zero local sketches because geometry is imported. Selected fiducial/channel voids, tube lumen, lid holes, arm clearance and service-view liquid terminations were visually reviewed. These establish executable CAD and selected visible features, not comprehensive interference, motion, tolerance, optical, structural, thermal or fluid validation.

## Software/document package
- `IMPLEMENTATION.md`: proposed acquisition/control state machine, identity/endpoint contract, fault response and scientific assumptions. Not deployed software.
- `REFERENCE_ANALYSIS.md`: unexecuted Python reference for finite-state stationary distributions, committors, MFPTs, fluxes and coarse-grained entropy, with analytic fixtures. Not original-engine integration.
- `VALIDATION.md`: proposed acceptance and test gates; no fabricated physical/software passes.
- `BOM.md`: CAD occurrences, not procurement release.
- `CHECKLIST.md`: achieved CAD scope versus remaining engineering.

## Remaining engineering / release blockers
1. Complete rated pressure/relief, gas/exhaust and inline sensor plumbing; liquid fittings/strain relief; electrical harnesses and hardware interlocks.
2. Detail second-axis stage drive, carriages, focus and travel stops. Freeze vendor-specific optics, interfaces and alignment mechanisms.
3. Complete structural/door/equipment retention, grounding and seals, then full interference/tolerance and relevant physical analyses.
4. Select component ratings and approved proof/leak procedures. **No operating pressure is assigned. A controller's possible7 bar capability is not a chip or assembly rating.**
5. Validate chip geometry, fabrication, wetting and viability for the selected cell type and assay.
6. Implement/run control, tracking, RNA processing, rate estimation, uncertainty analysis and original Waddington adapter against real code/API and data. No experimental result has been generated.
7. Obtain applicable independent safety, containment, ethics and regulatory assessment. No compliance or performance claim is made.
