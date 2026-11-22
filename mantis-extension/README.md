# Thermodynamic Waddington for Mantis

A Mantis panel for maps built from a Thermodynamic Waddington export
(`thermodynamic_waddington.export.export_cells`, `tw-export`, or
`benchmarks/export_examples.py` in the pipeline repository). Each point is a cell; its `tw_*`
fields carry the committor, commitment phase, free energy, and calibrated irreversibility.

The panel reads those fields from the active map and shows:

- **Where the lineage commits**: the first stage whose mean committor exceeds 0.5.
- **Landscape, three views**, all computed in the panel from the map's own cells:
  - *Cells*: every cell at its `tw_x` / `tw_y` position with the commitment-flow web behind
    it, coloured by committor, free energy, irreversibility or (when the export carries
    `tw_committor_sd`, from `tw-export --ensemble-fits N`) uncertainty, the committor's spread
    across independent refits; start, commitment and end stages
    labelled. Hovering a cell on the map rings it; opening a cell draws its path to commitment.
  - *Surface 3D*: the Waddington surface, a density pseudo-potential U = -ln(p/p_max + 0.04)
    coloured by the committor, with the q = 0.5 separatrix and the commitment stage marked.
    Drag to turn, scroll to zoom, double-click to reset.
  - *Commitment field*: the committor as a field (kernel regression), streamlines of its
    gradient, isocommittors every 0.2 and the q = 0.5 separatrix.
- **Readouts**: commitment stage, transition state, irreversibility and lineage length.
- **Commitment barrier**: the potential of mean force along the committor,
  F(q) = -kT ln P(q), computed in the panel with the pipeline's own definition (Gaussian KDE,
  60-point grid, maximum over the 5-95% interior), so it matches the Python result on the same
  cells, with a 95% band from resampling the cells and a chart of the barrier against the number
  of cells. It rises with cell number (the bandwidth shrinks as n^-1/5), so the panel also gives
  the resolution-matched barrier (fixed h = 0.08 in committor units, reflected at the
  boundaries, bootstrap bias-corrected), which matches the export's `barrier_resolved_kt` and
  does not drift with cell number. Single-fit barriers understate their uncertainty about
  twofold; quote them with the spread across independent fits.
- **Commitment flow between stages**: downstream edges counted between annotated stages, the
  share running forward along the committor, and the largest stage-to-stage flows.
- **Irreversibility by phase**: the per-cell noise-matched excess as box plots per phase, with
  the share of cells in each tail; balanced tails are what noise alone gives.
- **Commitment phases** and a **stage × phase** heatmap: the share of each stage's cells in
  each phase, with each stage's noise-matched irreversibility (about zero means reversible
  within the data's own noise).
- **The open cell**: its place on the committor, free energy, irreversibility, its path to
  commitment (from each cell, the downstream neighbour furthest along the committor), and its
  commitment flow drawn on the map as edges to the neighbours further along the committor.

A cell opens when it is clicked on the map or selected on its own.

**Motion.** The landscape is live: in the cells view, pulses travel from each cell toward the
neighbours further along its committor; on the 3D surface they stream downhill as it turns; in
the field view, particles follow the committor's gradient. On opening, the barrier curve draws
itself in and the phases, heatmap and flows fill in. The **Live** button pauses and resumes the
motion; it is off when the system asks for reduced motion. Nothing animates while the landscape is
scrolled out of view or the tab is hidden.

Actions write to the map: **Create phase bags** makes one bag per commitment phase
(`TW: progenitor` ... `TW: committed`), **Highlight transition** emphasises that bag,
**Select transition** selects its cells, and **Clear flow** removes the overlay. Closing the
panel releases the overlay and the emphasis.

The panel is near-black with the gradient of its headline as light behind frosted cards, and
adapts to its column: one column when narrow, four readouts across from 600 px, and a
two-column layout from 640 px. Drag the panel's splitter in Mantis to give it room.

## Package

```text
tw-mantis-extension/
  mantis.extension.json
  dist/
    panel.js
    panel.css
```

## Permissions

| Permission | Needed for |
| --- | --- |
| `maps:read` | `maps.getActive`, `maps.getPoints`, and the `point.clicked` / `point.hovered` events |
| `selection:read` | the `selection.changed` event (a cell selected on its own opens in the panel) |
| `bags:write` | `bags.create` (one bag per commitment phase) |
| `selection:write` | `bags.highlight`, `bags.unhighlight`, `selection.set` |
| `edges:write` | `edges.set` and `edges.clear` (the commitment-flow overlay) |

No host script, backend, commands or network access.

## Data

Create a space from `<prefix>_cells.csv`, or from an `.h5ad` whose `obs` carries the `tw_*`
columns and whose `obsm['X_umap']` holds the coordinates. Every `tw_*` column must be kept as a
field so it reaches each point's `metadata`. The panel only reads; it never changes the fields.

## Install

```bash
cd tw-mantis-extension
zip -r thermodynamic-waddington.mantisx mantis.extension.json dist README.md
```

In Mantis: open the space, **Verticals** menu, **Import Extension**, choose the `.mantisx`,
review the trust dialog, install. The panel appears under **Extensions** as
**Thermodynamic Waddington**.
