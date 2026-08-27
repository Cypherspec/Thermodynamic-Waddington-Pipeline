from __future__ import annotations

import html
import json
import math
from pathlib import Path
from typing import Any, Sequence

from .arrays import mean, percentile
from .cell_atlas import interpret_population
from .model import LandscapeFit
from .multimodal import cell_modalities

SAND = "#f5f5f7"
PAPER = "#ffffff"
INK = "#1d1d1f"
MUTED = "#6e6e73"
LINE = "#d2d2d7"
TEAL = "#0071e3"
GOLD = "#86868b"


def _finite(values: Sequence[float]) -> list[float]:
    return [float(value) for value in values if math.isfinite(float(value))]


def _safe_range(values: Sequence[float]) -> tuple[float, float]:
    finite = _finite(values)
    if not finite:
        return 0.0, 1.0
    low, high = min(finite), max(finite)
    return low, high if high > low else low + 1.0


def _norm(value: float, low: float, high: float) -> float:
    return max(0.0, min(1.0, (float(value) - low) / max(1e-12, high - low)))


def _rgb(color: str) -> tuple[int, int, int]:
    return tuple(int(color[index:index + 2], 16) for index in (1, 3, 5))


def _energy_color(value: float, low: float, high: float) -> str:
    ratio = _norm(value, low, high)
    stops = [(0.0, "#245c59"), (0.22, "#4c8178"), (0.45, "#a9b59a"), (0.62, "#dec993"), (0.82, "#c87f61"), (1.0, "#91483f")]
    for (left, left_color), (right, right_color) in zip(stops, stops[1:]):
        if ratio <= right:
            amount = (ratio - left) / max(1e-12, right - left)
            a, b = _rgb(left_color), _rgb(right_color)
            return "#%02x%02x%02x" % tuple(round(x + (y - x) * amount) for x, y in zip(a, b))
    return stops[-1][1]


def _embedding_bounds(fit: LandscapeFit) -> tuple[float, float, float, float]:
    points = fit.embedding or [[0.0, 0.0]]
    xs = [float(point[0]) for point in points]
    ys = [float(point[1]) for point in points]
    x0, x1 = _safe_range(xs)
    y0, y1 = _safe_range(ys)
    return x0 - (x1 - x0) * 0.055, x1 + (x1 - x0) * 0.055, y0 - (y1 - y0) * 0.055, y1 + (y1 - y0) * 0.055


def _grid_values(fit: LandscapeFit, cols: int = 24, rows: int = 18) -> tuple[list[list[float]], tuple[float, float, float, float]]:
    points = fit.embedding or [[0.0, 0.0]]
    energies = [float(value) for value in fit.energies] or [0.0]
    bounds = _embedding_bounds(fit)
    x0, x1, y0, y1 = bounds
    observed = [(float(point[0]), float(point[1]), energy) for point, energy in zip(points, energies)]
    distances: list[float] = []
    for index, (x, y, _) in enumerate(observed):
        nearest = sorted(math.hypot(x - ox, y - oy) for other, (ox, oy, _) in enumerate(observed) if other != index)
        distances.extend(nearest[:3])
    bandwidth = max(1e-9, percentile(distances, 0.55) * 1.8 if distances else 1.0)
    grid: list[list[float]] = []
    for row in range(rows + 1):
        y = y0 + (y1 - y0) * row / rows
        output_row: list[float] = []
        for col in range(cols + 1):
            x = x0 + (x1 - x0) * col / cols
            weighted = [(math.exp(-((px - x) ** 2 + (py - y) ** 2) / (2 * bandwidth * bandwidth)), energy) for px, py, energy in observed]
            denominator = sum(weight for weight, _ in weighted)
            output_row.append(sum(weight * energy for weight, energy in weighted) / denominator if denominator else mean(energies))
        grid.append(output_row)
    return grid, bounds


def _contour_segments(grid: Sequence[Sequence[float]], level: float, x0: float, x1: float, y0: float, y1: float) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    rows, cols = len(grid) - 1, len(grid[0]) - 1
    result: list[tuple[tuple[float, float], tuple[float, float]]] = []
    for row in range(rows):
        for col in range(cols):
            corners = [(col / cols, row / rows, grid[row][col]), ((col + 1) / cols, row / rows, grid[row][col + 1]), ((col + 1) / cols, (row + 1) / rows, grid[row + 1][col + 1]), (col / cols, (row + 1) / rows, grid[row + 1][col])]
            hits: list[tuple[float, float]] = []
            for first, second in ((corners[0], corners[1]), (corners[1], corners[2]), (corners[2], corners[3]), (corners[3], corners[0])):
                a, b = first[2], second[2]
                if (a < level and b < level) or (a > level and b > level) or a == b:
                    continue
                fraction = (level - a) / (b - a)
                hits.append((first[0] + fraction * (second[0] - first[0]), first[1] + fraction * (second[1] - first[1])))
            for first, second in zip(hits[::2], hits[1::2]):
                result.append(((x0 + first[0] * (x1 - x0), y0 + first[1] * (y1 - y0)), (x0 + second[0] * (x1 - x0), y0 + second[1] * (y1 - y0))))
    return result


def _basin_map(fit: LandscapeFit) -> dict[int, dict[str, Any]]:
    result: dict[int, dict[str, Any]] = {}
    records = fit.diagnostics.get("basins", [])
    for basin_index, record in enumerate(records if isinstance(records, list) else []):
        if not isinstance(record, dict):
            continue
        for member in record.get("members", []):
            result[int(member)] = {"basin": basin_index, "mass": float(record.get("mass", 0.0)), "minimum_energy": float(record.get("minimum_energy", 0.0))}
    return result


def _edge_stats(fit: LandscapeFit) -> dict[int, dict[str, Any]]:
    result: dict[int, dict[str, Any]] = {}
    for edge in fit.edges:
        source = int(edge["source"])
        target = int(edge["target"])
        item = result.setdefault(source, {"outgoing": 0, "alignment": 0.0, "barrier": 0.0, "targets": []})
        item["outgoing"] += 1
        item["alignment"] += float(edge.get("alignment", 0.0))
        item["barrier"] = max(float(item["barrier"]), float(edge.get("work", 0.0)))
        item["targets"].append(target)
    for item in result.values():
        item["alignment"] /= max(1, item["outgoing"])
    return result


def _phenotype_payload(fit: LandscapeFit, index: int) -> dict[str, Any] | None:
    for item in fit.diagnostics.get("cell_phenotypes", []) if isinstance(fit.diagnostics.get("cell_phenotypes", []), list) else []:
        if isinstance(item, dict) and int(item.get("cell_index", -1)) == index:
            return item
    return None


def _atlas_payload(fit: LandscapeFit, index: int) -> dict[str, Any] | None:
    expression = fit.metadata.get("expression", [])
    if not isinstance(expression, list) or index >= len(expression):
        return None
    interpretations = fit.diagnostics.get("cell_interpretations")
    if isinstance(interpretations, list) and index < len(interpretations) and isinstance(interpretations[index], dict):
        return interpretations[index]
    return None


def _cell_payload(fit: LandscapeFit):
    energies = [float(value) for value in fit.energies]
    uncertainties = [float(value) for value in (fit.uncertainties or [0.0] * len(energies))]
    energy_low, energy_high = _safe_range(energies)
    uncertainty_low, uncertainty_high = _safe_range(uncertainties)
    basins, edges = _basin_map(fit), _edge_stats(fit)
    expression, velocity, names = fit.metadata.get("expression", []), fit.metadata.get("velocity", []), fit.metadata.get("gene_names", [])
    payload: list[dict[str, Any]] = []
    for index, (point, energy) in enumerate(zip(fit.embedding, energies)):
        uncertainty = uncertainties[index] if index < len(uncertainties) else 0.0
        edge = edges.get(index, {"outgoing": 0, "alignment": 0.0, "barrier": 0.0, "targets": []})
        basin = basins.get(index, {})
        position, uncertainty_position = _norm(energy, energy_low, energy_high), _norm(uncertainty, uncertainty_low, uncertainty_high)
        state_class = "deep basin / low effective energy" if position <= 0.24 else "ridge / activated state" if position >= 0.76 else "transition corridor / metastable state"
        confidence = "uncertain local estimate" if uncertainty_position >= 0.72 else "high-confidence local estimate" if uncertainty_position <= 0.28 else "moderate-confidence local estimate"
        flow = float(edge["alignment"])
        direction = "strongly aligned with local downhill flow" if flow >= 0.55 else "weakly aligned with local flow" if flow >= 0.2 else "opposes the inferred local flow" if flow <= -0.2 else "locally transverse to the inferred flow"
        payload.append({
            "index": index,
            "id": str(fit.cell_ids[index]) if index < len(fit.cell_ids) else f"cell_{index:05d}",
            "label": str(fit.labels[index]) if index < len(fit.labels) else "unlabeled",
            "gene_names": [str(item) for item in names] if isinstance(names, list) else [],
            "gene_values": list(expression[index]) if isinstance(expression, list) and index < len(expression) and isinstance(expression[index], list) else [],
            "velocity_values": list(velocity[index]) if isinstance(velocity, list) and index < len(velocity) and isinstance(velocity[index], list) else [],
            "velocity_observed": bool(fit.diagnostics.get("velocity_observed", fit.metadata.get("velocity_observed", True))),
            "x": float(point[0]), "y": float(point[1]), "energy": energy, "energy_position": position,
            "uncertainty": uncertainty, "uncertainty_position": uncertainty_position,
            "outgoing_edges": int(edge["outgoing"]), "alignment": flow, "barrier_proxy": float(edge["barrier"]), "targets": edge["targets"],
            "basin": basin.get("basin"), "basin_mass": basin.get("mass"), "basin_minimum": basin.get("minimum_energy"),
            "state_class": state_class, "confidence": confidence, "direction": direction,
            "phenotype": _phenotype_payload(fit, index), "atlas": _atlas_payload(fit, index),
            "scenes": [scene for scene in fit.diagnostics.get("live_scenes", []) if isinstance(scene, dict) and int(scene.get("cell_index", -1)) == index],
            "modalities": cell_modalities(fit.metadata.get("modalities", {}), index),
            "measurement_ledger": (fit.diagnostics.get("cell_measurement_ledgers", [])[index] if isinstance(fit.diagnostics.get("cell_measurement_ledgers"), list) and index < len(fit.diagnostics.get("cell_measurement_ledgers", [])) else []),
            "population_context": (fit.diagnostics.get("cell_population_context", [])[index] if isinstance(fit.diagnostics.get("cell_population_context"), list) and index < len(fit.diagnostics.get("cell_population_context", [])) else []),
        })
    return payload


def _project(x: float, y: float, bounds: tuple[float, float, float, float], width: float, height: float) -> tuple[float, float]:
    x0, x1, y0, y1 = bounds
    return 58 + _norm(x, x0, x1) * (width - 116), height - 58 - _norm(y, y0, y1) * (height - 116)


def _surface_svg(fit: LandscapeFit, width: int = 1100, height: int = 700) -> str:
    grid, bounds = _grid_values(fit)
    low, high = percentile([value for row in grid for value in row], 0.02), percentile([value for row in grid for value in row], 0.98)
    rows, cols = len(grid) - 1, len(grid[0]) - 1
    polygons: list[str] = []
    for row in range(rows):
        for col in range(cols):
            values = [grid[row][col], grid[row][col + 1], grid[row + 1][col + 1], grid[row + 1][col]]
            points = []
            for px, py in ((col, row), (col + 1, row), (col + 1, row + 1), (col, row + 1)):
                x0, x1, y0, y1 = bounds
                points.append(_project(x0 + px / cols * (x1 - x0), y0 + py / rows * (y1 - y0), bounds, width, height))
            polygons.append(f'<polygon points="{" ".join(f"{x:.1f},{y:.1f}" for x,y in points)}" fill="{_energy_color(mean(values), low, high)}" fill-opacity=".94" stroke="#fffaf2" stroke-opacity=".16" stroke-width=".7"><title>F / kT = {mean(values):.3f}</title></polygon>')
    contours: list[str] = []
    for fraction in (.12, .24, .36, .48, .60, .72, .84):
        level = low + fraction * (high - low)
        for first, second in _contour_segments(grid, level, *bounds):
            a, b = _project(*first, bounds, width, height), _project(*second, bounds, width, height)
            contours.append(f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{b[0]:.1f}" y2="{b[1]:.1f}" class="contour"/>')
    cells = _cell_payload(fit)
    by_index = {cell["index"]: cell for cell in cells}
    flow: list[str] = []
    for edge in fit.edges[::max(1, len(fit.edges) // 220)]:
        source, target = int(edge["source"]), int(edge["target"])
        if source not in by_index or target not in by_index:
            continue
        a, b = _project(by_index[source]["x"], by_index[source]["y"], bounds, width, height), _project(by_index[target]["x"], by_index[target]["y"], bounds, width, height)
        flow.append(f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{b[0]:.1f}" y2="{b[1]:.1f}" class="flow"/>')
    nodes: list[str] = []
    for cell in cells:
        x, y = _project(cell["x"], cell["y"], bounds, width, height)
        nodes.append(f'<circle class="cell" data-index="{cell["index"]}" cx="{x:.1f}" cy="{y:.1f}" r="{4 + cell["uncertainty_position"] * 4:.1f}" fill="{_energy_color(cell["energy"], low, high)}"><title>{html.escape(cell["id"])} · F={cell["energy"]:.3f} kT</title></circle>')
    return f'<svg id="landscape-map" viewBox="0 0 {width} {height}" role="img" aria-label="2D thermodynamic RNA free-energy topographical landscape"><defs><linearGradient id="ground2d" x2="0" y2="1"><stop stop-color="#fbfaf7"/><stop offset="1" stop-color="#e8dfd2"/></linearGradient></defs><rect x="20" y="20" width="1060" height="650" rx="26" fill="url(#ground2d)"/><g id="map2d-world"><g class="terrain">{"".join(polygons)}</g><g class="contours">{"".join(contours)}</g><g class="flow">{"".join(flow)}</g><g class="cells">{"".join(nodes)}</g></g><text x="58" y="53" class="map-title">F(x) / kT · 2D topographical RNA landscape</text><text x="58" y="75" class="map-subtitle">basin depth · contour ridges · directed velocity flow · click a cell</text><text x="550" y="674" class="axis">embedding coordinate 1</text><text x="30" y="380" class="axis" transform="rotate(-90 30 380)">embedding coordinate 2</text></svg>'


def _height_mesh_svg(fit: LandscapeFit, width: int = 1100, height: int = 620) -> str:
    grid, bounds = _grid_values(fit)
    flat = [value for row in grid for value in row]
    low, high = percentile(flat, .02), percentile(flat, .98)
    terrain = {"grid": grid, "bounds": bounds, "low": low, "high": high, "width": width, "height": height}
    encoded = html.escape(json.dumps(terrain, separators=(",", ":")), quote=True)
    return f'<div class="terrain3d-wrap"><canvas id="landscape-3d" width="{width}" height="{height}" data-terrain="{encoded}" role="img" aria-label="3D thermodynamic RNA free-energy terrain"></canvas><div class="terrain3d-legend">drag to orbit · scroll to zoom · click a cell · double-click to reset · vertical axis is effective F / kT</div></div>'


def _insights(fit: LandscapeFit, cells: Sequence[dict[str, Any]]) -> list[str]:
    if not cells:
        return ["No cells were available."]
    deepest, highest, uncertain = min(cells, key=lambda item: item["energy"]), max(cells, key=lambda item: item["energy"]), max(cells, key=lambda item: item["uncertainty"])
    basin_counts: dict[int, int] = {}
    for cell in cells:
        if cell.get("basin") is not None:
            basin_counts[int(cell["basin"])] = basin_counts.get(int(cell["basin"]), 0) + 1
    largest = max(basin_counts, key=basin_counts.get) if basin_counts else None
    alignment_value = fit.diagnostics.get("velocity_alignment", fit.diagnostics.get("mean_velocity_alignment"))
    alignment = float(alignment_value) if alignment_value is not None else 0.0
    return [f"Deepest sampled cell: {deepest['id']} at F={deepest['energy']:.3f} kT; this is a sampled basin point, not proof of commitment.", f"Highest sampled ridge: {highest['id']} at F={highest['energy']:.3f} kT.", f"Largest uncertainty: {uncertain['id']} with σ={uncertain['uncertainty']:.3f}.", f"Mean directed RNA-velocity alignment is {alignment:.3f}; arrows are directional evidence, not proof of conservative dynamics." if fit.diagnostics.get("velocity_observed", True) else "RNA velocity was not observed; directed-flow alignment is unavailable and not interpreted as measured dynamics.", f"The largest inferred basin contains {basin_counts[largest]} cells." if largest is not None else "No populated basins were assigned."]


def _metric_card(title: str, value: object, detail: str) -> str:
    return f'<div class="metric"><div class="metric-title">{html.escape(title)}</div><div class="metric-value">{html.escape(str(value))}</div><div class="metric-detail">{html.escape(detail)}</div></div>'


def render_report(fit: LandscapeFit, output: str | Path) -> None:
    output = Path(output)
    cells = _cell_payload(fit)
    diagnostics = fit.diagnostics
    coverage = float(diagnostics.get("coverage", 0.0))
    velocity_observed = bool(diagnostics.get("velocity_observed", fit.metadata.get("velocity_observed", True)))
    velocity_status = str(diagnostics.get("velocity_status", fit.metadata.get("velocity_status", "observed")))
    alignment_value = diagnostics.get("velocity_alignment", diagnostics.get("mean_velocity_alignment"))
    alignment = float(alignment_value) if alignment_value is not None else float("nan")
    rows = "".join(f'<tr><td>{html.escape(str(label))}</td><td>{mean(values):+.3f}</td><td>{len(values)}</td></tr>' for label, values in sorted(((str(label), [float(energy) for energy, item in zip(fit.energies, fit.labels) if str(item) == str(label)]) for label in set(fit.labels)), key=lambda item: mean(item[1])))
    compact = {key: diagnostics[key] for key in ("coverage", "velocity_status", "velocity_observed", "velocity_alignment", "reversibility_gap", "advanced_metrics", "thermodynamic_summary", "geometric_summary", "warnings", "phenotype_summary") if key in diagnostics}
    data_json = json.dumps(cells, separators=(",", ":"))
    script = r'''function esc(value){return String(value).replace(/[&<>"']/g,function(c){return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c];});}
function fmt(value){var n=Number(value);return Number.isFinite(n)?n.toFixed(3):'n/a';}
function velocityText(cell){return cell.velocity_observed===false?'not observed':fmt(cell.alignment);}
function featureTable(cell){var names=cell.gene_names||[],values=cell.gene_values||[],velocity=cell.velocity_values||[];if(!names.length&&values.length)names=values.map(function(_,i){return 'feature_'+i;});var tiles=names.map(function(name,i){var value=Number(values[i]||0);var intensity=Math.max(0,Math.min(255,Math.round(128+value*42)));var color=value>=0?'rgb(190,'+(150-intensity*.18)+','+(115-intensity*.12)+')':'rgb('+(70-intensity*.12)+','+(120-intensity*.08)+','+(150+intensity*.1)+')';return '<div class="heat-tile" title="'+esc(name)+' · expression '+fmt(value)+'" style="background:'+color+'"><span>'+esc(name)+'</span><b>'+fmt(value)+'</b></div>';}).join('');return '<div class="gene-table"><div class="readout-title">complete expression and velocity data · '+names.length+' features</div><div class="cell-heatmap">'+tiles+'</div><table><thead><tr><th>gene / feature</th><th>RNA expression</th><th>RNA velocity</th></tr></thead><tbody>'+names.map(function(name,i){return '<tr><td>'+esc(name)+'</td><td>'+fmt(values[i])+'</td><td>'+fmt(velocity[i])+'</td></tr>';}).join('')+'</tbody></table></div>';}
function programRows(items){return (items||[]).map(function(item){var raw=item.score===null?item.strength:item.score;var matched=item.matched_markers||item.basis||[];var missing=item.missing_markers||[];var z=Number(item.z_score);var available=Number.isFinite(Number(raw));var zAvailable=Number.isFinite(z);var evidence=item.evidence||'inferred';var detail=(item.output||item.statement||item.description||'RNA-inferred program')+' · '+evidence+' · '+matched.length+' matched / '+missing.length+' missing markers';var shown=available?fmt(raw):(zAvailable?fmt(z)+' z':(matched.length?'RNA signal present':'not identifiable from supplied genes'));return '<div class="program-row"><div><b>'+esc(item.name||item.domain||item.key||'program')+'</b><small>'+esc(detail)+'</small><small>'+(matched.length?'matched: '+esc(matched.slice(0,8).join(', ')):'no marker names matched')+(missing.length?' · missing: '+esc(missing.slice(0,8).join(', ')):'')+'</small></div><strong class="'+(available||zAvailable?'':'na')+'">'+shown+'</strong></div>';}).join('');}
function modalityMarkup(m){if(!m)return '';var names=Object.keys(m);return '<div class="atlas-block"><div class="readout-title">multimodal evidence ledger</div><div class="atlas-caution">RNA remains the thermodynamic Waddington field. Other modalities are displayed only when supplied as measured data; missing modalities are not inferred.</div>'+names.map(function(name){var item=m[name]||{},values=item.values,features=item.feature_names||[];var measured=Array.isArray(values);var preview=measured?values.slice(0,12).map(function(v,i){return esc((features[i]||'feature_'+i)+': '+fmt(v));}).join(' · '):'not supplied';return '<div class="program-row"><div><b>'+esc(name.replaceAll('_',' '))+'</b><small>'+esc(item.source||'unknown source')+' · '+(measured?'measured':'unavailable')+'</small></div><strong class="'+(measured?'':'na')+'">'+preview+'</strong></div>';}).join('')+'</div>';}
function ledgerMarkup(ledger,context){var rows=(ledger||[]).map(function(item){var top=(item.top_features||[]).map(function(feature){return esc(feature.feature)+' '+fmt(feature.value);}).join(' · ')||'no finite values';return '<div class="program-row"><div><b>'+esc(String(item.modality||'').replaceAll('_',' '))+'</b><small>'+esc(item.source||'')+' · '+esc(item.caveat||'')+'</small></div><strong class="'+(item.status==='measured'?'':'na')+'">'+esc(item.status||'unavailable')+' · '+top+'</strong></div>';}).join('');var deviations=(context||[]).map(function(item){var outliers=(item.outliers||[]).slice(0,4).map(function(feature){return esc(feature.feature)+' z='+fmt(feature.z);}).join(' · ')||'no ranked deviations';return '<div class="program-row"><div><b>'+esc(String(item.modality||'').replaceAll('_',' '))+' · population context</b><small>'+esc(item.interpretation||'')+'</small></div><strong>'+fmt(item.population_distance)+' · '+outliers+'</strong></div>';}).join('');if(!rows&&!deviations)return '';return '<div class="atlas-block"><div class="readout-title">measured modality ledger and population context</div><div class="atlas-caution">Measured modality values are shown exactly when supplied. Population context is a descriptive comparison, not a clinical or mechanistic conclusion.</div>'+rows+deviations+'</div>';}
function phenotypeMarkup(p){if(!p)return '<div class="atlas-caution">No phenotype layer was produced.</div>';var outputs=(p.dominant_outputs||[]);if(!outputs.length&&p.evidence)outputs=p.evidence.filter(function(x){return x.domain==='function';}).map(function(x){return x.statement;});return '<div class="atlas-block"><div class="readout-title">predicted cell state and function</div><div class="readout-grid"><div class="readout-stat"><span>candidate identity</span><b>'+esc(p.identity||'unresolved')+' · '+fmt(p.identity_confidence)+'</b></div><div class="readout-stat"><span>developmental stage</span><b>'+esc(p.developmental_stage||'unresolved')+' · '+fmt(p.stage_confidence)+'</b></div><div class="readout-stat"><span>damage / stress</span><b>'+esc(p.damage_state||'unresolved')+' · '+fmt(p.damage_score)+'</b></div><div class="readout-stat"><span>viability evidence</span><b>'+esc(p.viability_state||'unresolved')+' · '+fmt(p.viability_score)+'</b></div><div class="readout-stat"><span>cell cycle</span><b>'+esc(p.cycle_state||'unresolved')+'</b></div><div class="readout-stat"><span>dominant outputs</span><b>'+esc(outputs.slice(0,4).join(' · ')||'unavailable')+'</b></div></div>'+programRows(p.evidence)+'<div class="atlas-caution">'+esc((p.caveats||[])[0]||'RNA-program inference; validate independently.')+'</div></div>';}
function atlasMarkup(a){if(!a)return '';var caveat=(a.cautions||[]).join(' · ')||'RNA-program inference; validate with orthogonal measurements.';return '<div class="atlas-block"><div class="readout-title">RNA-derived organelle and functional program map</div><div class="atlas-caution">'+esc(a.evidence_level||'program-level inference')+' · measured: '+esc((a.measured_modalities||[]).join(', '))+' · inferred: '+esc((a.inferred_modalities||[]).join(', '))+'</div>'+programRows(a.scores)+'<div class="atlas-caution">'+esc(caveat)+'</div></div>';}
function proxyMap(cell){var names=['mitochondrion','nucleus','ER','Golgi','lysosome','ribosome','membrane'];var scores=(cell.atlas&&cell.atlas.scores)||[];var values=names.map(function(name){var item=scores.find(function(x){return String(x.compartment||'').toLowerCase().indexOf(name.toLowerCase())>=0||String(x.name||'').toLowerCase().indexOf(name.toLowerCase())>=0;});return item&&Number.isFinite(Number(item.z_score))?Number(item.z_score):0;});var max=Math.max(1,...values.map(function(v){return Math.abs(v);}));return '<div class="proxy-map"><div class="readout-title">individual cell RNA program topology · evidence-scaled</div><div class="proxy-cell"><div class="proxy-core">RNA<br>state</div>'+names.map(function(name,i){var angle=i*360/names.length;var strength=Math.min(1,Math.abs(values[i])/max);return '<div class="proxy-organelle" style="--angle:'+angle+'deg;--level:'+strength+';--signal:'+(values[i]>=0?'#326b67':'#9b5a50')+'"><span>'+esc(name)+'</span><b>'+fmt(values[i])+' z</b></div>';}).join('')+'</div><div class="atlas-caution">Each spoke is an RNA marker-program z-score. Positive and negative signal are shown; zero means no identifiable signal in the supplied gene set, not absence of an organelle. This is not microscopy.</div></div>';}
function sceneMarkup(cell){var scenes=cell.scenes||[];if(!scenes.length)return '';return '<div class="atlas-block"><div class="readout-title">live-mode evidence scenes</div>'+scenes.map(function(scene){return '<div class="program-row"><div><b>scene '+esc(scene.number)+': '+esc(scene.title)+'</b><small>'+esc(scene.narration)+'</small><small>'+esc(scene.grounding||'Computed from supplied inputs')+'</small></div><strong>'+esc((scene.evidence||[]).join(' · '))+'</strong></div>';}).join('')+'</div>';}
function selectCell(index){var cell=CELLS.find(function(item){return item.index===Number(index);});if(!cell)return;window.__selectedCell=cell.index;document.querySelectorAll('.cell').forEach(function(node){node.classList.toggle('selected',Number(node.dataset.index)===cell.index);});select.value=String(cell.index);document.getElementById('selected-cell-badge').textContent=cell.id+' · '+cell.label;document.getElementById('readout').scrollIntoView({behavior:'smooth',block:'start'});readout.classList.remove('readout-enter');void readout.offsetWidth;readout.classList.add('readout-enter');readout.innerHTML='<div class="readout-title">cell interpretation · complete multi-layer atlas</div><h3>'+esc(cell.id)+'</h3><div class="cell-subtitle">'+esc(cell.label)+' · '+esc(cell.state_class)+'</div><div class="readout-grid"><div class="readout-stat"><span>effective energy F / kT</span><b>'+fmt(cell.energy)+'</b></div><div class="readout-stat"><span>bootstrap uncertainty</span><b>σ '+fmt(cell.uncertainty)+'</b></div><div class="readout-stat"><span>basin / attractor</span><b>'+esc(cell.basin===null?'unassigned':'basin '+cell.basin)+(cell.attractor?' · candidate well':'')+'</b></div><div class="readout-stat"><span>transition pressure</span><b>'+fmt(cell.barrier_proxy)+'</b></div><div class="readout-stat"><span>flow alignment</span><b>'+velocityText(cell)+'</b></div><div class="readout-stat"><span>confidence</span><b>'+esc(cell.confidence)+'</b></div></div><p class="interpretation">'+esc(cell.state_class)+'. The local RNA velocity is '+esc(cell.velocity_observed===false?'not observed in this dataset; flow-dependent interpretation is unavailable':cell.direction)+'. This remains the original thermodynamic RNA Waddington landscape: F(x)/kT is an effective inference from expression and velocity, not a molecular equilibrium state function.</p>'+phenotypeMarkup(cell.phenotype)+atlasMarkup(cell.atlas)+proxyMap(cell)+modalityMarkup(cell.modalities)+ledgerMarkup(cell.measurement_ledger,cell.population_context)+sceneMarkup(cell)+featureTable(cell);document.getElementById('complete-cell-data').innerHTML=featureTable(cell);if(window.__draw3d)window.__draw3d();}
function appendMapPoint(cell){var svg=document.getElementById('landscape-map');if(!svg)return;var b=window.__mapBounds||[-1,1,-1,1],x0=58+(Number(cell.x)-b[0])/Math.max(1e-9,b[1]-b[0])*984,y0=642-(Number(cell.y)-b[2])/Math.max(1e-9,b[3]-b[2])*584,c=document.createElementNS('http://www.w3.org/2000/svg','circle');c.setAttribute('class','cell');c.dataset.index=cell.index;c.setAttribute('cx',x0);c.setAttribute('cy',y0);c.setAttribute('r','6');c.setAttribute('fill','#c69b57');c.setAttribute('tabindex','0');c.setAttribute('aria-label',cell.id);c.addEventListener('click',function(){selectCell(cell.index);});c.addEventListener('mouseenter',function(){var hover=document.getElementById('hover-readout');if(hover)hover.textContent=cell.id+' · F / kT '+fmt(cell.energy)+' · '+cell.state_class;});c.addEventListener('keydown',function(e){if(e.key==='Enter'||e.key===' ')selectCell(cell.index);});svg.querySelector('.cells').appendChild(c);}
function addUploadedRows(payload,filename){var matrix=payload.expression||payload.X||[];var names=payload.gene_names||payload.genes||[];if(!Array.isArray(matrix)||!matrix.length){alert('Expected JSON with expression: [[...],[...]].');return;}var all=matrix.map(function(row){return row.map(Number).filter(Number.isFinite);});if(!names.length)names=all[0].map(function(_,i){return 'feature_'+i;});var maxAbs=Math.max(1,...all.flat().map(Math.abs));var start=CELLS.length;var rows=all.map(function(values,i){var half=Math.max(1,Math.floor(values.length/2));var x=values.slice(0,half).reduce(function(a,b){return a+b;},0)/Math.max(1,half)/maxAbs;var y=values.slice(half).reduce(function(a,b){return a+b;},0)/Math.max(1,values.length-half)/maxAbs;var energy=values.reduce(function(a,b){return a+Math.abs(b);},0)/Math.max(1,values.length);return {index:start+i,id:(payload.cell_ids||[])[i]||filename+'_'+i,label:(payload.labels||[])[i]||'uploaded RNA',gene_names:names,gene_values:values,velocity_values:(payload.velocity||[])[i]||[],velocity_observed:Array.isArray(payload.velocity),x:Math.max(-1,Math.min(1,x)),y:Math.max(-1,Math.min(1,y)),energy:energy,uncertainty:0,basin:null,attractor:false,barrier_proxy:0,alignment:0,confidence:'local upload preview',state_class:'uploaded RNA state · refit for authoritative landscape',direction:'not yet integrated',outgoing_edges:0,phenotype:null,atlas:null};});CELLS.push.apply(CELLS,rows);rows.forEach(function(cell){appendMapPoint(cell);var option=document.createElement('option');option.value=cell.index;option.textContent=cell.id+' · uploaded RNA preview';select.appendChild(option);});if(window.__draw3d)window.__draw3d();document.getElementById('upload-status').textContent='Added '+rows.length+' cell(s) to the interactive map. Preview heatmaps are live; refit the bundle to compute thermodynamic energy, basins, flow, and bootstrap uncertainty.';}
function install2D(){var svg=document.getElementById('landscape-map'),layer=svg&&svg.querySelector('#map2d-world');if(!svg||!layer)return;var state={zoom:1,x:0,y:0,drag:false,lastX:0,lastY:0};function clampPan(){var limit=Math.max(0,(state.zoom-1)*470);state.x=Math.max(-limit,Math.min(limit,state.x));state.y=Math.max(-limit*.66,Math.min(limit*.66,state.y));}function apply(){clampPan();layer.style.transform='translate('+state.x+'px,'+state.y+'px) scale('+state.zoom+')';svg.dataset.zoom=state.zoom.toFixed(2);}function reset(){state={zoom:1,x:0,y:0,drag:false,lastX:0,lastY:0};apply();}svg.addEventListener('wheel',function(e){e.preventDefault();var box=svg.getBoundingClientRect(),mx=e.clientX-box.left,my=e.clientY-box.top,old=state.zoom,next=Math.max(1,Math.min(4,old*(e.deltaY<0?1.12:.89)));state.x=state.x-(mx-box.width/2)*(next/old-1);state.y=state.y-(my-box.height/2)*(next/old-1);state.zoom=next;apply();},{passive:false});svg.addEventListener('pointerdown',function(e){if(state.zoom<=1)return;if(e.target.closest&&e.target.closest('.cell'))return;state.drag=true;state.lastX=e.clientX;state.lastY=e.clientY;svg.setPointerCapture(e.pointerId);});svg.addEventListener('pointermove',function(e){if(!state.drag)return;state.x+=e.clientX-state.lastX;state.y+=e.clientY-state.lastY;state.lastX=e.clientX;state.lastY=e.clientY;apply();});['pointerup','pointercancel','pointerleave'].forEach(function(name){svg.addEventListener(name,function(){state.drag=false;});});svg.addEventListener('dblclick',reset);window.__reset2d=reset;apply();}
function install3D(){var canvas=document.getElementById('landscape-3d');if(!canvas)return;var ctx=canvas.getContext('2d'),terrain=JSON.parse(canvas.dataset.terrain||'{}'),state={zoom:1,x:0,y:0,rx:-.62,ry:-.62,drag:false,lastX:0,lastY:0};function project(ix,iy,z){var grid=terrain.grid||[],rows=grid.length,cols=rows?(grid[0]||[]).length:0,low=Number(terrain.low||0),high=Number(terrain.high||1),xx=(ix/Math.max(1,cols-1)-.5)*760,yy=(iy/Math.max(1,rows-1)-.5)*410,zz=(Number(z)-low)/Math.max(1e-9,high-low)*290,c=Math.cos(state.ry),s=Math.sin(state.ry),xr=xx*c-yy*s,yr=xx*s+yy*c,pitch=Math.cos(state.rx),sp=Math.sin(state.rx);return [canvas.width/2+(xr+state.x)*state.zoom,canvas.height*.63+(yr*pitch-zz*sp+state.y)*state.zoom];}function color(t){t=Math.max(0,Math.min(1,t));var a,b,u;if(t<.5){a=[31,91,107];b=[38,151,145];u=t*2;}else{a=[38,151,145];b=[232,174,82];u=(t-.5)*2;}return 'rgb('+Math.round(a[0]+(b[0]-a[0])*u)+','+Math.round(a[1]+(b[1]-a[1])*u)+','+Math.round(a[2]+(b[2]-a[2])*u)+')';}function draw(){var w=canvas.width,h=canvas.height,g=terrain.grid||[],rows=g.length,cols=rows?(g[0]||[]).length:0,low=Number(terrain.low||0),high=Number(terrain.high||1);ctx.clearRect(0,0,w,h);var bg=ctx.createLinearGradient(0,0,0,h);bg.addColorStop(0,'#fbfbfd');bg.addColorStop(1,'#e8e1d7');ctx.fillStyle=bg;ctx.fillRect(0,0,w,h);ctx.lineJoin='round';for(var r=rows-2;r>=0;r--){for(var c=0;c<cols-1;c++){var z00=g[r][c],z10=g[r][c+1],z11=g[r+1][c+1],z01=g[r+1][c],points=[project(c,r,z00),project(c+1,r,z10),project(c+1,r+1,z11),project(c,r+1,z01)],v=(z00+z10+z11+z01)/4,t=(v-low)/Math.max(1e-9,high-low);ctx.beginPath();ctx.moveTo(points[0][0],points[0][1]);points.slice(1).forEach(function(q){ctx.lineTo(q[0],q[1]);});ctx.closePath();ctx.fillStyle=color(t);ctx.fill();ctx.strokeStyle='rgba(255,255,255,.28)';ctx.lineWidth=.7;ctx.stroke();}}ctx.strokeStyle='#4b5563';ctx.lineWidth=1.5;ctx.beginPath();ctx.moveTo(68,h-66);ctx.lineTo(w-72,h-66);ctx.moveTo(68,h-66);ctx.lineTo(218,h-250);ctx.moveTo(68,h-66);ctx.lineTo(68,42);ctx.stroke();ctx.fillStyle='#27313a';ctx.font='600 italic 20px Georgia';ctx.fillText('effective RNA free-energy terrain',48,30);ctx.font='12px ui-monospace';ctx.fillText('embedding x / y',w-190,h-38);ctx.fillText('F / kT',24,85);(CELLS||[]).forEach(function(cell){if(cell.hidden)return;var b=terrain.bounds||[-1,1,-1,1],rows2=Math.max(1,rows-1),cols2=Math.max(1,cols-1),ix=(cell.x-b[0])/Math.max(1e-9,b[1]-b[0])*cols2,iy=(cell.y-b[2])/Math.max(1e-9,b[3]-b[2])*rows2,q=project(ix,iy,cell.energy),active=Number(cell.index)===Number(window.__selectedCell);ctx.beginPath();ctx.arc(q[0],q[1],active?9:5,0,Math.PI*2);ctx.fillStyle=active?'#f05a5a':'#f4d06f';ctx.fill();ctx.strokeStyle='#ffffff';ctx.lineWidth=active?3:1;ctx.stroke();});}function reset(){state={zoom:1,x:0,y:0,rx:-.62,ry:-.62,drag:false,lastX:0,lastY:0};draw();}canvas.addEventListener('wheel',function(e){e.preventDefault();state.zoom=Math.max(1,Math.min(4,state.zoom*(e.deltaY<0?1.12:.89)));draw();},{passive:false});canvas.addEventListener('pointerdown',function(e){state.drag=true;state.lastX=e.clientX;state.lastY=e.clientY;canvas.setPointerCapture(e.pointerId);});canvas.addEventListener('pointermove',function(e){if(!state.drag)return;var dx=e.clientX-state.lastX,dy=e.clientY-state.lastY;state.ry+=dx*.009;state.rx=Math.max(-1.45,Math.min(.25,state.rx-dy*.009));state.lastX=e.clientX;state.lastY=e.clientY;draw();});['pointerup','pointercancel','pointerleave'].forEach(function(name){canvas.addEventListener(name,function(){state.drag=false;});});canvas.addEventListener('dblclick',reset);canvas.addEventListener('click',function(e){var rect=canvas.getBoundingClientRect(),scaleX=canvas.width/Math.max(1,rect.width),scaleY=canvas.height/Math.max(1,rect.height),px=(e.clientX-rect.left)*scaleX,py=(e.clientY-rect.top)*scaleY,best=null,dist=Infinity;(CELLS||[]).forEach(function(cell){if(cell.hidden)return;var b=terrain.bounds||[-1,1,-1,1],rows2=Math.max(1,(terrain.grid||[]).length-1),cols2=Math.max(1,(terrain.grid||[[]])[0].length-1),ix=(cell.x-b[0])/Math.max(1e-9,b[1]-b[0])*cols2,iy=(cell.y-b[2])/Math.max(1e-9,b[3]-b[2])*rows2,q=project(ix,iy,cell.energy),d=(q[0]-px)*(q[0]-px)+(q[1]-py)*(q[1]-py);if(d<dist){dist=d;best=cell;}});if(best&&dist<900)selectCell(best.index);});window.__draw3d=draw;window.__reset3d=reset;draw();}
function installCellSearch(){var input=document.getElementById('cell-search'),clear=document.getElementById('clear-search');if(!input)return;function apply(){var query=String(input.value||'').trim().toLowerCase();CELLS.forEach(function(cell){var hay=[cell.id,cell.label,cell.state_class,cell.confidence,cell.basin===null?'':('basin '+cell.basin)].concat(cell.gene_names||[]).join(' ').toLowerCase();cell.hidden=!!query&&hay.indexOf(query)<0;var node=document.querySelector('.cell[data-index="'+cell.index+'"]');if(node)node.style.display=cell.hidden?'none':'';});if(window.__draw3d)window.__draw3d();document.getElementById('search-status').textContent=query?'showing '+CELLS.filter(function(cell){return !cell.hidden;}).length+' of '+CELLS.length+' cells':'showing all '+CELLS.length+' cells';}input.addEventListener('input',apply);if(clear)clear.addEventListener('click',function(){input.value='';apply();input.focus();});apply();}
function addUploadHandlers(){var input=document.getElementById('bundle-upload');if(!input)return;input.addEventListener('change',function(){Array.prototype.forEach.call(input.files,function(file){var reader=new FileReader();reader.onload=function(){try{addUploadedRows(JSON.parse(reader.result),file.name);}catch(e){document.getElementById('upload-status').textContent='Could not read '+file.name+': '+e.message;}};reader.readAsText(file);});input.value='';});}
var CELLS=__DATA__,select=document.getElementById('cell-select'),readout=document.getElementById('readout'),mode=document.getElementById('map-mode');window.__mapBounds=__BOUNDS__;function setMode(value){mode.value=value;document.getElementById('map-shell-2d').hidden=value!=='2d';document.getElementById('map-shell-3d').hidden=value!=='3d';document.getElementById('map-shell-3d').classList.toggle('active-3d',value==='3d');if(value==='3d'&&window.__draw3d)requestAnimationFrame(window.__draw3d);}setMode('2d');mode.addEventListener('change',function(){setMode(mode.value);});CELLS.forEach(function(cell){var option=document.createElement('option');option.value=cell.index;option.textContent=cell.id+' · '+cell.label;select.appendChild(option);});document.querySelectorAll('.cell').forEach(function(node){node.addEventListener('click',function(){selectCell(node.dataset.index);});node.addEventListener('mouseenter',function(){var cell=CELLS.find(function(item){return item.index===Number(node.dataset.index);});var hover=document.getElementById('hover-readout');if(cell&&hover)hover.textContent=cell.id+' · F / kT '+fmt(cell.energy)+' · '+cell.state_class;});node.addEventListener('focus',function(){node.dispatchEvent(new Event('mouseenter'));});});select.addEventListener('change',function(){if(select.value!=='')selectCell(select.value);});document.addEventListener('click',function(event){var target=event.target.closest&&event.target.closest('.cell');if(target)selectCell(target.dataset.index);});document.getElementById('download-cell').addEventListener('click',function(){var cell=CELLS.find(function(item){return item.index===Number(select.value);});if(!cell)return;var link=document.createElement('a');link.href=URL.createObjectURL(new Blob([JSON.stringify(cell,null,2)],{type:'application/json'}));link.download=(cell.id||'cell')+'_thermodynamic_waddington.json';link.click();});document.getElementById('reset').addEventListener('click',function(){document.querySelectorAll('.cell').forEach(function(node){node.classList.remove('selected');});select.value='';window.__selectedCell=null;readout.innerHTML='<div class="placeholder">Click a point on the landscape or choose a cell above.</div>';document.getElementById('complete-cell-data').innerHTML='<div class="atlas-caution">Select a cell to see the individual RNA heatmap and complete feature table.</div>';});addUploadHandlers();install2D();install3D();installCellSearch();
var liveIndex=0,liveTime=0,livePlaying=false,liveFrame=0,liveLast=0,liveSceneIndex=0;
function liveCell(){return CELLS[liveIndex%Math.max(1,CELLS.length)];}
function renderLive(){if(!CELLS.length)return;var cell=liveCell(),scene=(cell.scenes||[])[liveSceneIndex%(Math.max(1,(cell.scenes||[]).length))];document.getElementById('live-title').textContent='RNA field replay · '+cell.id;document.getElementById('live-detail').textContent=(scene?scene.title+' · '+scene.narration:'F / kT '+fmt(cell.energy)+' · '+cell.state_class)+' · basin '+(cell.basin===null?'unassigned':cell.basin)+' · '+(cell.velocity_observed===false?'velocity not observed':'flow alignment '+fmt(cell.alignment));document.getElementById('live-time').textContent='t = '+liveTime.toFixed(2)+' · cell '+(liveIndex%CELLS.length+1)+' / '+CELLS.length;document.getElementById('live-scrub').value=String(Math.round((liveIndex%CELLS.length)/Math.max(1,CELLS.length-1)*100));document.querySelector('.live-panel').style.setProperty('--live-progress',String((liveIndex%CELLS.length)/Math.max(1,CELLS.length-1)));document.querySelectorAll('.cell').forEach(function(node){node.classList.toggle('live-focus',Number(node.dataset.index)===Number(cell.index));});window.__liveMode=true;window.__selectedCell=cell.index;if(window.__draw3d)window.__draw3d();var evidence=document.getElementById('live-evidence');if(evidence){var direction=cell.velocity_observed===false?'velocity unavailable':'RNA velocity '+velocityText(cell);evidence.innerHTML='<span class="live-chip">scene '+(liveIndex+1)+' / '+CELLS.length+'</span><span class="live-chip">'+esc(cell.id)+'</span><span class="live-chip">F '+fmt(cell.energy)+' kT</span><span class="live-chip">'+esc(direction)+'</span><span class="live-chip">basin '+esc(cell.basin===null?'unassigned':cell.basin)+'</span>'+(scene?'<span class="live-chip">'+esc(scene.title)+'</span>':'');}}
function liveRenderSelected(){var cell=liveCell();if(cell){selectCell(cell.index);}renderLive();}
function liveStep(){if(!CELLS.length)return;liveIndex=(liveIndex+1)%CELLS.length;liveSceneIndex=0;liveTime+=.25;liveRenderSelected();}
function liveLoop(now){if(!livePlaying)return;if(!liveLast)liveLast=now;var elapsed=now-liveLast,speed=Number(document.getElementById('live-speed').value||1),interval=Math.max(180,900/speed);if(elapsed>=interval){liveLast=now;liveStep();}liveFrame=requestAnimationFrame(liveLoop);}
document.getElementById('live-toggle').addEventListener('click',function(){livePlaying=!livePlaying;this.textContent=livePlaying?'pause live mode':'play live mode';liveLast=0;if(livePlaying)liveFrame=requestAnimationFrame(liveLoop);else cancelAnimationFrame(liveFrame);});document.getElementById('live-step').addEventListener('click',liveStep);document.getElementById('live-scrub').addEventListener('input',function(){liveIndex=Math.round(Number(this.value)/100*Math.max(0,CELLS.length-1));liveTime=liveIndex*.25;liveRenderSelected();});renderLive();
'''.replace('__DATA__', data_json).replace('__BOUNDS__', json.dumps(list(_embedding_bounds(fit)), separators=(',', ':')))
    template = '''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Thermodynamic Waddington · Cell Fate Terrain</title><style>:root{--sand:__SAND__;--paper:__PAPER__;--ink:__INK__;--muted:__MUTED__;--line:__LINE__;--teal:__TEAL__;--gold:__GOLD__;font-family:-apple-system,BlinkMacSystemFont,"SF Pro Display","Helvetica Neue",Arial,sans-serif;color:var(--ink);background:var(--sand);font-synthesis:none}*{box-sizing:border-box}body{margin:0;background:linear-gradient(180deg,#fbfbfd 0%,#f5f5f7 100%);min-height:100vh;-webkit-font-smoothing:antialiased}main{max-width:1280px;margin:auto;padding:52px 28px 96px}.eyebrow{font:600 10px ui-monospace,SFMono-Regular,Menlo,monospace;letter-spacing:.16em;color:var(--teal);text-transform:uppercase}h1,h2,h3{font-family:-apple-system,BlinkMacSystemFont,"SF Pro Display","Helvetica Neue",Arial,sans-serif;font-style:normal;font-weight:600;letter-spacing:-.035em}h1{font-size:clamp(48px,8vw,92px);line-height:.94;margin:18px 0 22px}h2{font-size:28px;margin:28px 0 15px}.dek,.interpretation,.cell-subtitle{font-family:-apple-system,BlinkMacSystemFont,"SF Pro Text","Helvetica Neue",Arial,sans-serif;color:#515154;line-height:1.55}.dek{font-size:19px;max-width:780px}.notice{border:1px solid var(--line);border-radius:16px;padding:16px 18px;margin:18px 0;color:#515154;background:#ffffffcc;font-size:13px;line-height:1.6;box-shadow:0 8px 24px #1d1d1f0a}.panel{background:var(--paper);border:1px solid var(--line);border-radius:24px;box-shadow:0 14px 40px #1d1d1f0d;padding:16px;overflow:hidden}.surface{padding:0}.map-panel{padding:0;overflow:hidden}.map-shell{position:relative;width:100%;overflow:hidden;border-radius:18px;background:#f2f2f7}.map-shell[hidden]{display:none}.map-shell svg{display:block;width:100%;height:auto;touch-action:none;cursor:grab;transform-origin:center center}.map-shell svg:active{cursor:grabbing}.map-shell svg text{user-select:none}.map-shell canvas{display:block;width:100%;height:auto;border-radius:18px;background:#f2f2f7}.map-shell .contour{stroke:#ffffff;stroke-opacity:.68;stroke-width:1.1;fill:none}.map-shell .flow{stroke:#1d1d1f;stroke-opacity:.18;stroke-width:1;marker-end:url(#arrow)}.map-shell .cell{cursor:pointer;stroke:#ffffff;stroke-width:1;transition:stroke-width .18s,filter .18s,opacity .18s}.map-shell .cell:hover,.map-shell .mesh-cell:hover,.map-shell .selected{stroke:#0071e3;stroke-width:3;filter:drop-shadow(0 0 5px #0071e366)}.map-shell .mesh-stem{stroke:#8e8e93;stroke-width:1;stroke-opacity:.55}.live-panel{margin-top:14px;position:relative;background:linear-gradient(135deg,#ffffff,#f2f2f7);overflow:hidden}.live-panel:before{content:'';display:block;height:3px;width:calc(var(--live-progress,0)*100%);background:linear-gradient(90deg,#0071e3,#8e8e93);transition:width .45s ease}.live-stage{display:flex;align-items:center;gap:18px;min-height:84px}.live-stage p{margin:6px 0 0;color:var(--muted);font-size:12px}.live-orb{width:54px;height:54px;border-radius:50%;background:radial-gradient(circle at 35% 30%,#ffffff,#8e8e93 48%,#0071e3 82%);box-shadow:0 0 0 7px #ffffffcc,0 0 34px #0071e333;animation:live-pulse 2.1s ease-in-out infinite}@keyframes live-pulse{50%{transform:scale(1.08);box-shadow:0 0 0 13px #ffffff99,0 0 45px #0071e355}}.zoom-hint{font:10px ui-monospace,SFMono-Regular,Menlo,monospace;color:var(--muted);margin:7px 2px}.selected-badge{font:11px ui-monospace,SFMono-Regular,Menlo,monospace;color:var(--teal);border:1px solid var(--line);background:#ffffff;padding:9px 12px;border-radius:999px}.controls{display:flex;flex-wrap:wrap;gap:8px;margin:13px 0;align-items:center}button,.select,.upload-control{background:#ffffff;border:1px solid var(--line);color:var(--ink);padding:10px 14px;border-radius:999px;font:500 11px ui-monospace,SFMono-Regular,Menlo,monospace;cursor:pointer;transition:background .18s,border-color .18s,transform .18s,box-shadow .18s}button:hover,.upload-control:hover{border-color:var(--teal);box-shadow:0 5px 16px #0071e31f}button:active{transform:scale(.98)}.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(165px,1fr));gap:10px;margin:20px 0}.metric{background:#f5f5f7;border:1px solid var(--line);border-radius:16px;padding:15px}.metric-title{font:600 10px ui-monospace,SFMono-Regular,Menlo,monospace;color:var(--muted);text-transform:uppercase;letter-spacing:.09em}.metric-value{font-size:25px;margin:7px 0}.metric-detail{font-size:12px;color:var(--muted)}.grid{display:grid;grid-template-columns:minmax(0,1.8fr) minmax(280px,.8fr);gap:18px;align-items:start}.readout{min-height:500px}.placeholder{color:var(--muted);font-family:-apple-system,BlinkMacSystemFont,"SF Pro Text","Helvetica Neue",Arial,sans-serif;padding:23px}.readout-title{font:600 10px ui-monospace,SFMono-Regular,Menlo,monospace;letter-spacing:.12em;color:var(--teal);text-transform:uppercase;margin-bottom:10px}.readout-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(145px,1fr));gap:8px;margin:15px 0}.readout-stat{border:1px solid var(--line);border-radius:14px;padding:11px;background:#ffffff}.readout-stat span{display:block;font:600 10px ui-monospace,SFMono-Regular,Menlo,monospace;color:var(--muted);text-transform:uppercase}.readout-stat b{display:block;margin-top:6px;font-size:13px}.program-row{display:flex;justify-content:space-between;gap:12px;border-top:1px solid var(--line);padding:10px 0}.program-row small{display:block;color:var(--muted);margin-top:4px;max-width:660px}.program-row strong{white-space:nowrap}.program-row .na{color:var(--muted);font-weight:400}.atlas-block{margin-top:19px;padding-top:15px;border-top:1px solid var(--line)}.atlas-caution{font-size:12px;color:#515154;background:#f5f5f7;border-left:3px solid var(--teal);padding:10px 12px;margin-top:12px;line-height:1.5}.cell-heatmap{display:grid;grid-template-columns:repeat(auto-fill,minmax(110px,1fr));gap:5px;padding:9px;background:#f5f5f7;border:1px solid var(--line);max-height:300px;overflow:auto}.heat-tile{min-height:48px;border-radius:8px;padding:5px;color:#ffffff;text-shadow:0 1px 2px #1d1d1f88;font:10px ui-monospace,SFMono-Regular,Menlo,monospace;display:flex;flex-direction:column;justify-content:space-between}.heat-tile b{font-size:12px}.proxy-map{margin-top:18px;padding-top:15px;border-top:1px solid var(--line)}.proxy-cell{width:250px;height:210px;margin:10px auto;position:relative;border:2px solid #0071e3;border-radius:50%;background:radial-gradient(circle at 50% 50%,#d2d2d7 0 17%,#f5f5f7 18% 32%,#d1e3f5 33% 62%,#a7c7e6 63% 78%,#5b8fc4 79%)}.proxy-core{position:absolute;left:92px;top:76px;width:66px;height:48px;border-radius:50%;background:#6e6e73;color:#ffffff;text-align:center;padding-top:12px;font:10px ui-monospace,SFMono-Regular,Menlo,monospace}.proxy-organelle{position:absolute;left:calc(50% - 42px);top:calc(50% - 17px);width:84px;height:34px;transform:rotate(var(--angle)) translateY(-84px) rotate(calc(var(--angle) * -1));font:9px ui-monospace,SFMono-Regular,Menlo,monospace;text-align:center;color:#1d1d1f}.proxy-organelle span,.proxy-organelle b{display:block;background:#ffffffdd;border-radius:7px;padding:2px}.proxy-organelle b{background:#0071e3;color:#ffffff}.gene-table{max-height:560px;overflow:auto;border:1px solid var(--line);margin-top:18px;border-radius:14px}.gene-table table{border-collapse:collapse;width:100%;font:11px ui-monospace,SFMono-Regular,Menlo,monospace}.gene-table th,.gene-table td{text-align:left;padding:6px 8px;border-bottom:1px solid #e5e5ea}.gene-table th{position:sticky;top:0;background:#f5f5f7}.side-list{margin:0;padding-left:19px;line-height:1.75;color:#515154}.table-wrap{overflow:auto}.label-table{width:100%;border-collapse:collapse;font:12px ui-monospace,SFMono-Regular,Menlo,monospace}.label-table th,.label-table td{text-align:left;padding:8px;border-bottom:1px solid var(--line)}.footnote{font:10px ui-monospace,SFMono-Regular,Menlo,monospace;color:var(--muted);margin-top:22px}@media(max-width:850px){.grid{grid-template-columns:1fr}h1{font-size:56px}}
.live-evidence{display:flex;gap:8px;flex-wrap:wrap;margin-top:12px;min-height:28px}.live-chip{display:inline-flex;align-items:center;border:1px solid var(--line);background:#ffffffcc;border-radius:999px;padding:7px 10px;font:10px ui-monospace,SFMono-Regular,Menlo,monospace;color:#515154;transition:transform .35s ease,opacity .35s ease}.live-focus{animation:cell-focus 1.2s ease-in-out infinite;transform-origin:center}@keyframes cell-focus{50%{filter:drop-shadow(0 0 11px #0071e399);stroke-width:4}}}</style></head><body><main><div class="eyebrow">THERMODYNAMIC WADDINGTON · EFFECTIVE RNA LANDSCAPE</div><h1>Cell fate, as terrain.</h1><p class="dek">The original thermodynamic RNA landscape: an effective F(x)/kT field from expression and directed RNA velocity, with basin structure, uncertainty, flow, and cell-state program interpretation.</p><div class="notice"><b>Data provenance.</b> {velocity_status}{" · measured RNA velocity is available for directed flow." if velocity_observed else " · no measured RNA velocity was supplied; velocity-dependent results are unavailable or provisional."}</div><div class="notice"><b>Interpretation boundary.</b> Damage, viability, organelle activity, identity, and stage are RNA-program inferences—not direct observations of organelles, proteins, metabolites, morphology, or physical damage. The terrain is an effective diagnostic, not a calibrated molecular equilibrium state function.</div><div class="metrics">__METRICS__</div><section class="panel map-panel"><div class="controls"><label>map view <select id="map-mode" class="select"><option value="2d">2D contour topograph</option><option value="3d">3D height terrain</option></select></label><span class="zoom-hint">Scroll = zoom · drag = pan/orbit · double-click = reset · minimum zoom 1×</span></div><div class="map-shell" id="map-shell-2d">__SURFACE__</div><div class="map-shell" id="map-shell-3d" hidden>__MESH__</div><div id="hover-readout" class="hover-readout" aria-live="polite">Hover a cell to preview its RNA-derived state.</div><div class="zoom-hint">2D preserves embedding geometry. 3D adds a rendered z-axis for effective energy height. Both are the same thermodynamic RNA landscape.</div><section class="live-panel"><div class="readout-title">cinematic RNA field replay · live mode</div><div class="live-stage"><div class="live-orb"></div><div><div class="readout-title">live field replay</div><b id="live-title">RNA field replay</b><p id="live-detail">Step through the fitted thermodynamic RNA landscape cell by cell.</p><p id="live-time">t = 0.00</p></div></div><div class="controls"><button id="live-toggle">play live mode</button><button id="live-step">step one cell</button><input id="live-scrub" type="range" min="0" max="100" value="0" style="flex:1;min-width:160px"><label class="zoom-hint">speed <input id="live-speed" type="range" min="0.25" max="3" step="0.25" value="1"></label></div><div id="live-evidence" class="live-evidence"></div></section></section><section class="controls"><span id="selected-cell-badge" class="selected-badge">no cell selected</span><select id="cell-select" class="select"><option value="">select a cell…</option></select><input id="cell-search" class="select" type="search" placeholder="search cells, states, genes…"><button id="clear-search">clear</button><span id="search-status" class="zoom-hint"></span><button id="download-cell">download selected cell JSON</button><label class="upload-control">add expression / velocity JSON<input id="bundle-upload" type="file" accept=".json,application/json" multiple hidden></label><button id="reset">reset selection</button><span id="upload-status" class="zoom-hint"></span></section><div class="grid"><section class="panel readout" id="readout"><div class="placeholder">Click a point on the landscape or choose a cell above. The individual cell heatmap and full RNA interpretation will appear here.</div></section><aside class="panel"><h2>Key inferences</h2><ul class="side-list">__INSIGHTS__</ul><h2>Cell fate label summary</h2><div class="table-wrap"><table class="label-table"><thead><tr><th>label</th><th>mean F / kT</th><th>n</th></tr></thead><tbody>__ROWS__</tbody></table></div></aside></div><section class="panel" id="complete-cell-data"><h2>Complete cell data</h2><div class="placeholder">Select a cell to see the individual RNA heatmap, organelle-program proxy map, and every expression/velocity value.</div></section><section class="panel"><h2>Diagnostics</h2><details><summary>JSON diagnostics</summary><pre>__DIAGNOSTICS__</pre></details></section><div class="footnote">Sand-white edition · one map view at a time · original thermodynamic RNA Waddington field · all biology panels are transcript-level inference.</div></main><script>__SCRIPT__</script></body></html>'''
    metrics = "".join([_metric_card("cells", len(fit.energies), "single-cell observations"), _metric_card("directed flow", len(fit.edges), "local transitions"), _metric_card("coverage", f"{coverage:.1%}", "reachable from reference"), _metric_card("alignment", f"{alignment:.3f}" if velocity_observed else "n/a", "measured RNA velocity" if velocity_observed else "no measured velocity; provisional"), _metric_card("candidate wells", len(fit.attractors), "low-energy attractors"), _metric_card("mean σ", f"{mean(fit.uncertainties):.3f}", "bootstrap uncertainty")])
    report = template.replace("__SAND__", SAND).replace("__PAPER__", PAPER).replace("__INK__", INK).replace("__MUTED__", MUTED).replace("__LINE__", LINE).replace("__TEAL__", TEAL).replace("__GOLD__", GOLD).replace("__METRICS__", metrics).replace("__SURFACE__", _surface_svg(fit)).replace("__MESH__", _height_mesh_svg(fit)).replace("__INSIGHTS__", "".join(f"<li>{html.escape(item)}</li>" for item in _insights(fit, cells))).replace("__ROWS__", rows or '<tr><td colspan="3">No labels supplied.</td></tr>').replace("__DIAGNOSTICS__", html.escape(json.dumps(compact, indent=2, sort_keys=True))).replace("__SCRIPT__", script)
    output.write_text(report)
