# CellFlux instrument (hardware)

CellFlux is the physical companion to the Thermodynamic Waddington pipeline: a
benchtop microfluidic instrument, designed to capture single-cell RNA (spliced
and unspliced, for RNA velocity) under controlled perturbations, so the readout
can be fed straight into the software's thermodynamic analysis.

## Status: a design, not a built instrument

Nothing here has been fabricated or tested. It is a hardware design:

- **`cellflux/`** - a Zoo (KCL 2.0) parametric CAD project, 91 parts: a droplet
  flow-focusing chip (real underside channels, ~30 um nozzle), the enclosure
  (base plate, uprights, roof, windowed panels, access doors), optics mounts
  (camera coupler and riser, filter cube), fluidics (supply lines, ferrules,
  flow sensor), and a control-electronics housing. See `cellflux/stuff/BOM.md`,
  `VALIDATION.md`, and `CHECKLIST.md`.
- **`DEVICE_DESIGN.md`** - the full instrument proposal: assays, component
  choices, microfabrication, the data pipeline, a validation plan, risk
  analysis, and a 2-3 year build timeline.

The BOM states plainly that quantities are CAD occurrences and that materials and
vendors are pending. Treat this as an engineering design under development.

## How it connects to the software

The instrument is meant to close the loop the software analyzes:

```
cells + perturbation
      -> CellFlux  (microfluidics + imaging + scRNA-seq, spliced/unspliced)
      -> RNA velocity
      -> thermodynamic_waddington  (entropy production, free-energy landscape,
                                    committor, commitment barrier)
      -> a thermodynamic readout of the perturbation
```

The preregistered commitment test in
[`../PREREGISTRATION_commitment.md`](../PREREGISTRATION_commitment.md) is exactly
the kind of experiment this instrument is designed to run: perturb cells at the
transition state versus the progenitor and measure the fate shift the software
predicts.

## Opening the CAD

The parts are KCL (kcl-lang, Zoo Design Studio). Open the project in Zoo Design
Studio, or use the KCL toolchain, to view and edit the parametric models. Each
`.kcl` file is a single part; `geometry.kcl` and `fluidRouting.kcl` are shared
construction libraries.
