/*
 * Thermodynamic Waddington - a Mantis panel for maps built from a Thermodynamic
 * Waddington export (one point per cell, `tw_*` fields in each point's metadata).
 *
 * What the panel does, in the order a reader of the code meets it:
 *
 *   1. Read the active map once and collect the cells that carry `tw_committor`.
 *   2. Summarise them: where the lineage commits, the commitment barrier from
 *      the potential of mean force along the committor (the pipeline's own
 *      definition, so the numbers agree), cells per commitment phase, and each
 *      stage's calibrated irreversibility.
 *   3. Draw the committor as a field over the cells, with the q = 0.5 separatrix
 *      -- the line the lineage crosses to commit.
 *   4. Turn results back onto the map: commitment phases become bags, downstream
 *      flow becomes edges, and a clicked cell's path to commitment is traced.
 *
 * It is one self-contained script, loaded inside a sandboxed iframe. Everything
 * is computed in the browser from the point metadata -- there is no backend.
 */
(function () {
    'use strict';

    /* ------------------------------------------------------------------ util */

    const $ = (id) => document.getElementById(id);
    const KT = 'kT';

    function num(value) {
        const n = typeof value === 'number' ? value : parseFloat(value);
        return Number.isFinite(n) ? n : null;
    }
    function mean(xs) {
        return xs.length ? xs.reduce((a, b) => a + b, 0) / xs.length : NaN;
    }
    function std(xs) {
        if (xs.length < 2) return 0;
        const m = mean(xs);
        return Math.sqrt(xs.reduce((a, b) => a + (b - m) * (b - m), 0) / (xs.length - 1));
    }
    function quantile(sorted, p) {
        if (!sorted.length) return NaN;
        const i = (sorted.length - 1) * p;
        const lo = Math.floor(i), hi = Math.ceil(i);
        return lo === hi ? sorted[lo] : sorted[lo] + (sorted[hi] - sorted[lo]) * (i - lo);
    }
    function fmt(n, d) {
        if (n === null || n === undefined || Number.isNaN(n)) return '--';
        return Number(n).toFixed(d === undefined ? 2 : d);
    }
    function escapeHtml(s) {
        return String(s).replace(/[&<>"']/g, (c) => (
            { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
    }

    /* ------------------------------------------------------------- phases */

    /* Commitment phases are bands of the committor q, the probability a cell's
     * lineage reaches its terminal fate. The separatrix is q = 0.5: below it a
     * cell is more likely to return, above it more likely to commit. */
    const PHASES = [
        { name: 'Uncommitted', lo: 0.0, hi: 0.15, color: '#6b7a90' },
        { name: 'Priming', lo: 0.15, hi: 0.5, color: '#4f8bff' },
        { name: 'Committing', lo: 0.5, hi: 0.85, color: '#34d3ee' },
        { name: 'Committed', lo: 0.85, hi: 1.01, color: '#8b7bff' },
    ];
    function phaseOf(q) {
        for (const p of PHASES) if (q >= p.lo && q < p.hi) return p.name;
        return PHASES[PHASES.length - 1].name;
    }
    function phaseByName(s) {
        const p = PHASES.find((p) => p.name.toLowerCase() === String(s).toLowerCase());
        return p ? p.name : null;
    }
    function phaseColor(name) {
        const p = PHASES.find((p) => p.name === name);
        return p ? p.color : '#6b7a90';
    }

    /** One cell from a map point, or null when the point carries no Thermodynamic Waddington fields. */
    function cellOf(point) {
        const m = point && point.metadata ? point.metadata : null;
        if (!m) return null;
        const q = num(m.tw_committor);
        if (q === null) return null;
        const named = m.tw_phase ? phaseByName(String(m.tw_phase)) : null;
        return {
            id: point.id,
            twId: m.tw_id !== undefined ? String(m.tw_id) : point.data_row,
            stage: m.tw_stage !== undefined ? String(m.tw_stage) : 'unlabelled',
            q: q,
            phase: named || phaseOf(q),
            energy: num(m.tw_energy),
            excess: num(m.tw_cyclic_excess),
            sd: num(m.tw_committor_sd),
            downstream: typeof m.tw_downstream === 'string' && m.tw_downstream.length > 0
                ? m.tw_downstream.split('|') : [],
            x: num(m.tw_x),
            y: num(m.tw_y),
        };
    }

    /* ----------------------------------------------------- potential of mean force */

    /* F(q) = -kT ln P(q) along the committor, with the pipeline's own definition:
     * a Gaussian KDE of the committor values on a 60-point grid, shifted to a zero
     * minimum, and the barrier the maximum of F over the interior 5-95% quantile
     * band (so a few outliers at q = 0 or 1 do not set the height). Silverman's
     * bandwidth is used, so the barrier rises slowly with the number of cells, as
     * in the export. */
    function silverman(samples) {
        const n = samples.length;
        if (n < 2) return 0.1;
        const s = std(samples);
        const sorted = samples.slice().sort((a, b) => a - b);
        const iqr = quantile(sorted, 0.75) - quantile(sorted, 0.25);
        const spread = iqr > 0 ? Math.min(s, iqr / 1.349) : s;
        return Math.max(1e-3, 0.9 * (spread || 0.1) * Math.pow(n, -1 / 5));
    }
    function pmf(qs) {
        const grid = [];
        for (let i = 0; i < 60; i += 1) grid.push(i / 59);
        const h = silverman(qs);
        const inv = 1 / (qs.length * h * Math.sqrt(2 * Math.PI));
        const density = grid.map((g) => {
            let acc = 0;
            for (let j = 0; j < qs.length; j += 1) {
                const u = (g - qs[j]) / h;
                acc += Math.exp(-0.5 * u * u);
            }
            return acc * inv;
        });
        const F = density.map((p) => -Math.log(p + 1e-9));
        const sorted = qs.slice().sort((a, b) => a - b);
        const lo = quantile(sorted, 0.05), hi = quantile(sorted, 0.95);
        let fMin = Infinity;
        grid.forEach((g, i) => { if (g >= lo && g <= hi && F[i] < fMin) fMin = F[i]; });
        if (!Number.isFinite(fMin)) fMin = Math.min.apply(null, F);
        const shifted = F.map((f) => f - fMin);
        let barrier = 0;
        grid.forEach((g, i) => { if (g >= lo && g <= hi && shifted[i] > barrier) barrier = shifted[i]; });
        return { grid, F: shifted, barrier, h, interior: [lo, hi] };
    }

    /* --------------------------------------------------------- committor field */

    /* A short Turbo-like ramp: cool (uncommitted) through to warm (committed). */
    function ramp(t) {
        t = Math.max(0, Math.min(1, t));
        const stops = [
            [0.00, [30, 40, 70]], [0.25, [40, 110, 200]], [0.5, [40, 200, 210]],
            [0.75, [120, 210, 150]], [1.0, [150, 130, 250]],
        ];
        for (let i = 1; i < stops.length; i += 1) {
            if (t <= stops[i][0]) {
                const [a, ca] = stops[i - 1], [b, cb] = stops[i];
                const f = (t - a) / (b - a || 1);
                return `rgb(${ca.map((v, k) => Math.round(v + (cb[k] - v) * f)).join(',')})`;
            }
        }
        return 'rgb(150,130,250)';
    }

    /* Kernel regression of the committor over an (x, y) grid, so the field is a
     * smooth surface rather than scattered dots, and the q = 0.5 separatrix can
     * be drawn where the regressed value crosses one half. */
    function drawField(canvas, cells) {
        const placed = cells.filter((c) => c.x !== null && c.y !== null);
        const dpr = Math.min(2, window.devicePixelRatio || 1);
        const W = canvas.clientWidth || 360, H = canvas.clientHeight || 200;
        canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr);
        const ctx = canvas.getContext('2d');
        ctx.scale(dpr, dpr);
        ctx.clearRect(0, 0, W, H);
        if (placed.length < 3) return;
        let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
        for (const c of placed) {
            minX = Math.min(minX, c.x); maxX = Math.max(maxX, c.x);
            minY = Math.min(minY, c.y); maxY = Math.max(maxY, c.y);
        }
        const pad = 10;
        const sx = (x) => pad + (x - minX) / (maxX - minX || 1) * (W - 2 * pad);
        const sy = (y) => H - pad - (y - minY) / (maxY - minY || 1) * (H - 2 * pad);
        const spanX = (maxX - minX) || 1, spanY = (maxY - minY) || 1;
        const bw = 0.12 * Math.sqrt(spanX * spanY);
        const reg = (px, py) => {
            let wsum = 0, qsum = 0;
            for (const c of placed) {
                const dx = c.x - px, dy = c.y - py;
                const w = Math.exp(-(dx * dx + dy * dy) / (2 * bw * bw));
                wsum += w; qsum += w * c.q;
            }
            return wsum > 1e-9 ? qsum / wsum : null;
        };
        // Soft field: a coarse grid of regressed committor, drawn as blocks.
        const step = 10;
        for (let px = 0; px < W; px += step) {
            for (let py = 0; py < H; py += step) {
                const gx = minX + (px - pad) / (W - 2 * pad) * spanX;
                const gy = minY + (H - pad - py) / (H - 2 * pad) * spanY;
                const q = reg(gx, gy);
                if (q === null) continue;
                ctx.fillStyle = ramp(q);
                ctx.globalAlpha = 0.5;
                ctx.fillRect(px, py, step, step);
            }
        }
        ctx.globalAlpha = 1;
        // The q = 0.5 separatrix: grid points where the regression crosses one half.
        ctx.fillStyle = 'rgba(255,255,255,0.92)';
        for (let px = 0; px < W; px += 3) {
            for (let py = 0; py < H; py += 3) {
                const gx = minX + (px - pad) / (W - 2 * pad) * spanX;
                const gy = minY + (H - pad - py) / (H - 2 * pad) * spanY;
                const q = reg(gx, gy);
                if (q !== null && Math.abs(q - 0.5) < 0.012) ctx.fillRect(px, py, 1.6, 1.6);
            }
        }
        // The cells themselves, so density reads.
        for (const c of placed) {
            ctx.beginPath();
            ctx.arc(sx(c.x), sy(c.y), 1.7, 0, Math.PI * 2);
            ctx.fillStyle = ramp(c.q);
            ctx.globalAlpha = 0.9;
            ctx.fill();
        }
        ctx.globalAlpha = 1;
    }

    /** F(q) curve with the barrier marked. */
    function drawPmf(canvas, p) {
        const dpr = Math.min(2, window.devicePixelRatio || 1);
        const W = canvas.clientWidth || 360, H = canvas.clientHeight || 120;
        canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr);
        const ctx = canvas.getContext('2d');
        ctx.scale(dpr, dpr);
        ctx.clearRect(0, 0, W, H);
        const pad = 6, fMax = Math.max.apply(null, p.F) || 1;
        const sx = (q) => pad + q * (W - 2 * pad);
        const sy = (f) => H - pad - (f / fMax) * (H - 2 * pad);
        // separatrix marker
        ctx.strokeStyle = 'rgba(255,255,255,0.18)';
        ctx.setLineDash([3, 3]); ctx.beginPath();
        ctx.moveTo(sx(0.5), pad); ctx.lineTo(sx(0.5), H - pad); ctx.stroke();
        ctx.setLineDash([]);
        // the curve
        const grad = ctx.createLinearGradient(0, 0, W, 0);
        grad.addColorStop(0, '#4f8bff'); grad.addColorStop(0.5, '#34d3ee'); grad.addColorStop(1, '#8b7bff');
        ctx.strokeStyle = grad; ctx.lineWidth = 2; ctx.lineJoin = 'round'; ctx.beginPath();
        p.grid.forEach((q, i) => { const X = sx(q), Y = sy(p.F[i]); i ? ctx.lineTo(X, Y) : ctx.moveTo(X, Y); });
        ctx.stroke();
        // barrier dot
        let bi = 0; p.F.forEach((f, i) => { if (f > p.F[bi]) bi = i; });
        ctx.fillStyle = '#eaf0ff';
        ctx.beginPath(); ctx.arc(sx(p.grid[bi]), sy(p.F[bi]), 3, 0, Math.PI * 2); ctx.fill();
    }

    /* --------------------------------------------------------- code-wave + shader */

    const CODEWAVE_GLYPHS = '░▒▓│┤┼┐└┴┬├─◇◆<>/\\|=+*·:;xX0O#%&';
    function makeCodewave(rows, cols) {
        rows = rows || 12; cols = cols || 150;
        let seed = 0x1a2b3c;
        const rnd = () => ((seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff);
        const lines = [];
        for (let r = 0; r < rows; r += 1) {
            let line = '';
            for (let c = 0; c < cols; c += 1) {
                line += rnd() < 0.2 ? ' ' : CODEWAVE_GLYPHS[Math.floor(rnd() * CODEWAVE_GLYPHS.length)];
            }
            lines.push(line);
        }
        return lines.join('\n');
    }

    const VERT = 'attribute vec2 p;void main(){gl_Position=vec4(p,0.0,1.0);}';
    const FRAG = [
        'precision highp float;',
        'uniform float u_time; uniform vec2 u_res;',
        'float hash(vec2 p){p=fract(p*vec2(123.34,456.21));p+=dot(p,p+45.32);return fract(p.x*p.y);}',
        'float noise(vec2 p){vec2 i=floor(p),f=fract(p);vec2 u=f*f*(3.0-2.0*f);',
        ' return mix(mix(hash(i),hash(i+vec2(1.0,0.0)),u.x),mix(hash(i+vec2(0.0,1.0)),hash(i+vec2(1.0,1.0)),u.x),u.y);}',
        'float fbm(vec2 p){float v=0.0,a=0.5;for(int i=0;i<5;i++){v+=a*noise(p);p*=2.0;a*=0.5;}return v;}',
        'void main(){',
        ' vec2 uv=gl_FragCoord.xy/u_res.xy;',
        ' vec2 p=uv*vec2(u_res.x/max(u_res.y,1.0),1.0)*2.3;',
        ' float t=u_time*0.045;',
        ' vec2 q=vec2(fbm(p+t),fbm(p+vec2(5.2,1.3)-t));',
        ' vec2 r=vec2(fbm(p+3.0*q+vec2(1.7,9.2)+0.12*t),fbm(p+3.0*q+vec2(8.3,2.8)-0.10*t));',
        ' float f=fbm(p+3.0*r);',
        ' vec3 col=mix(vec3(0.02,0.04,0.11),vec3(0.13,0.40,0.96),clamp(f*1.3,0.0,1.0));',
        ' col=mix(col,vec3(0.26,0.80,1.0),clamp(r.x*r.x*0.95,0.0,1.0));',
        ' col=mix(col,vec3(0.44,0.34,0.96),clamp(q.y*0.42,0.0,1.0));',
        ' col+=0.03*hash(uv+t);',
        ' gl_FragColor=vec4(col,1.0);',
        '}',
    ].join('\n');

    /* The WebGL liquid backdrop. Self-contained, capped resolution, animates only
     * when motion is allowed, and leaves the glyph field behind it untouched if
     * WebGL is unavailable. */
    function mountShader(canvas) {
        if (!canvas) return;
        let gl;
        try {
            gl = canvas.getContext('webgl', { alpha: false, antialias: false, powerPreference: 'low-power' })
                || canvas.getContext('experimental-webgl');
        } catch (e) { gl = null; }
        if (!gl) return;
        const reduce = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
        const compile = (type, src) => {
            const s = gl.createShader(type);
            gl.shaderSource(s, src); gl.compileShader(s);
            return gl.getShaderParameter(s, gl.COMPILE_STATUS) ? s : null;
        };
        const vs = compile(gl.VERTEX_SHADER, VERT), fs = compile(gl.FRAGMENT_SHADER, FRAG);
        const prog = gl.createProgram();
        if (!vs || !fs || !prog) return;
        gl.attachShader(prog, vs); gl.attachShader(prog, fs); gl.linkProgram(prog);
        if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) return;
        gl.useProgram(prog);
        const buf = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, buf);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
        const loc = gl.getAttribLocation(prog, 'p');
        gl.enableVertexAttribArray(loc);
        gl.vertexAttribPointer(loc, 2, gl.FLOAT, false, 0, 0);
        const uTime = gl.getUniformLocation(prog, 'u_time');
        const uRes = gl.getUniformLocation(prog, 'u_res');
        const start = performance.now();
        let raf = 0;
        const resize = () => {
            const dpr = Math.min(1.25, window.devicePixelRatio || 1);
            const w = Math.max(1, Math.round(canvas.clientWidth * dpr));
            const h = Math.max(1, Math.round(canvas.clientHeight * dpr));
            if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; }
            gl.viewport(0, 0, canvas.width, canvas.height);
        };
        const render = () => {
            resize();
            gl.uniform1f(uTime, (performance.now() - start) / 1000);
            gl.uniform2f(uRes, canvas.width, canvas.height);
            gl.drawArrays(gl.TRIANGLES, 0, 3);
            if (!reduce) raf = requestAnimationFrame(render);
        };
        render();
        window.addEventListener('resize', () => { if (reduce) render(); });
    }

    /* ------------------------------------------------------------------ template */

    function kpi(key, value, note, tone) {
        return '<div class="kpi' + (tone ? ' ' + tone : '') + '"><span class="k">' + escapeHtml(key)
            + '</span><span class="v">' + value + '</span><span class="n">' + escapeHtml(note) + '</span></div>';
    }

    const TEMPLATE = [
        '<div class="tw">',
        '  <header class="hero">',
        '    <canvas class="shader" id="twShader" aria-hidden="true"></canvas>',
        '    <div class="codewave" id="twCodewave" aria-hidden="true"></div>',
        '    <div class="head">',
        '      <span class="eyebrow"><i></i>THERMODYNAMIC WADDINGTON</span>',
        '      <span class="pill" id="twCount">--</span>',
        '    </div>',
        '    <h1 class="title">Commitment, <span class="key">barrier and irreversibility</span></h1>',
        '    <p class="sub">Where a lineage crosses the separatrix, the barrier it climbs in kT, and whether its flow circulates beyond its own noise.</p>',
        '  </header>',
        '  <div class="body" id="twBody">',
        '    <section class="kpis" id="twKpis"></section>',
        '    <section class="card">',
        '      <div class="card-h"><h2>Committor field</h2><button class="ghost" id="twBags">Make phase bags</button></div>',
        '      <canvas class="field" id="twField"></canvas>',
        '      <p class="caption" id="twFieldCap">The committor over the cells; the bright line is the q = 0.5 separatrix, where commitment tips.</p>',
        '    </section>',
        '    <section class="card">',
        '      <div class="card-h"><h2>Potential of mean force</h2><span class="aside" id="twPmfAside"></span></div>',
        '      <canvas class="pmf" id="twPmf"></canvas>',
        '      <p class="caption">F(q) = &minus;kT ln P(q) along the committor (Gaussian KDE, Silverman bandwidth), shifted to a zero minimum. The barrier is its height over the interior.</p>',
        '    </section>',
        '    <section class="card">',
        '      <div class="card-h"><h2>Commitment phases</h2><span class="aside" id="twPhaseAside"></span></div>',
        '      <div class="phases" id="twPhases"></div>',
        '    </section>',
        '    <section class="card">',
        '      <div class="card-h"><h2>Stage irreversibility</h2><button class="ghost" id="twFlow">Draw flow on map</button></div>',
        '      <div class="stages" id="twStages"></div>',
        '      <p class="caption" id="twIrrCap"></p>',
        '    </section>',
        '    <section class="card">',
        '      <div class="card-h"><h2>Cell</h2></div>',
        '      <div class="cell" id="twCell"><span class="hint-dot"></span>Click a cell on the map to place it on the committor and trace its path to commitment.</div>',
        '    </section>',
        '  </div>',
        '  <div class="toast" id="twToast" hidden></div>',
        '</div>',
    ].join('\n');

    /* ------------------------------------------------------------------ app */

    let CELLS = [];
    let MAP_ID = null;

    function toast(msg) {
        const el = $('twToast');
        if (!el) return;
        el.textContent = msg; el.hidden = false;
        clearTimeout(toast._t);
        toast._t = setTimeout(() => { el.hidden = true; }, 2600);
    }

    function renderSummary(cells) {
        const qs = cells.map((c) => c.q);
        const p = pmf(qs);
        // Commitment point: the stage whose mean committor sits closest to 0.5.
        const byStage = groupBy(cells, (c) => c.stage);
        let commitStage = '--', best = Infinity;
        for (const [stage, list] of byStage) {
            const d = Math.abs(mean(list.map((c) => c.q)) - 0.5);
            if (d < best) { best = d; commitStage = stage; }
        }
        const excesses = cells.map((c) => c.excess).filter((e) => e !== null);
        const meanExcess = excesses.length ? mean(excesses) : null;
        const committed = cells.filter((c) => c.q >= 0.5).length;

        $('twCount').textContent = cells.length.toLocaleString() + ' cells';
        $('twKpis').innerHTML = [
            kpi('Cells', cells.length.toLocaleString(), (byStage.size) + ' stages'),
            kpi('Commits at', escapeHtml(commitStage), 'crosses q = 0.5'),
            kpi('Barrier', fmt(p.barrier, 2) + ' <span class="unit">' + KT + '</span>',
                'over q ∈ [' + fmt(p.interior[0], 2) + ', ' + fmt(p.interior[1], 2) + ']'),
            kpi('Committed', Math.round(100 * committed / cells.length) + '<span class="unit">%</span>',
                'q ≥ 0.5', committed / cells.length >= 0.5 ? 'good' : ''),
        ].join('');

        drawField($('twField'), cells);
        drawPmf($('twPmf'), p);
        $('twPmfAside').textContent = 'Silverman h = ' + fmt(p.h, 3);

        // Phases
        const byPhase = groupBy(cells, (c) => c.phase);
        const total = cells.length;
        $('twPhaseAside').textContent = byPhase.size + ' phases';
        $('twPhases').innerHTML = PHASES.filter((ph) => byPhase.has(ph.name)).map((ph) => {
            const n = byPhase.get(ph.name).length;
            const pct = Math.round(100 * n / total);
            return '<div class="phase" data-phase="' + escapeHtml(ph.name) + '">'
                + '<span class="dot" style="background:' + ph.color + '"></span>'
                + '<span class="pname">' + escapeHtml(ph.name) + '</span>'
                + '<span class="bar"><span style="width:' + pct + '%;background:' + ph.color + '"></span></span>'
                + '<span class="pn">' + n + '</span></div>';
        }).join('');

        // Stage irreversibility
        const stageRows = Array.from(byStage.entries()).map(([stage, list]) => {
            const ex = list.map((c) => c.excess).filter((e) => e !== null);
            return { stage, n: list.length, excess: ex.length ? mean(ex) : null, q: mean(list.map((c) => c.q)) };
        }).sort((a, b) => a.q - b.q);
        const exVals = stageRows.map((r) => r.excess).filter((e) => e !== null);
        const exMax = exVals.length ? Math.max(0.01, Math.max.apply(null, exVals.map(Math.abs))) : 0.01;
        $('twStages').innerHTML = stageRows.map((r) => {
            const e = r.excess;
            const w = e === null ? 0 : Math.min(100, Math.abs(e) / exMax * 100);
            const sign = e !== null && e > 0 ? 'pos' : 'neg';
            return '<div class="stage"><span class="sname">' + escapeHtml(r.stage) + '</span>'
                + '<span class="sbar"><span class="' + sign + '" style="width:' + w + '%"></span></span>'
                + '<span class="sn">' + (e === null ? '--' : (e > 0 ? '+' : '') + fmt(e, 3)) + '</span></div>';
        }).join('');
        if (meanExcess === null) {
            $('twIrrCap').textContent = 'No tw_cyclic_excess field on this map.';
        } else {
            const verdict = meanExcess > 0.01 ? 'circulates beyond its noise'
                : meanExcess < -0.01 ? 'below the noise-matched null' : 'reversible within noise';
            $('twIrrCap').innerHTML = 'Mean cyclic excess ' + (meanExcess > 0 ? '+' : '') + fmt(meanExcess, 3)
                + ' &middot; the lineage ' + verdict + '. Positive is net circulation above a noise-matched null.';
        }
    }

    function groupBy(xs, key) {
        const m = new Map();
        for (const x of xs) {
            const k = key(x);
            if (!m.has(k)) m.set(k, []);
            m.get(k).push(x);
        }
        return m;
    }

    async function makePhaseBags() {
        const byPhase = groupBy(CELLS, (c) => c.phase);
        let made = 0;
        for (const ph of PHASES) {
            if (!byPhase.has(ph.name)) continue;
            const ids = byPhase.get(ph.name).map((c) => c.id);
            try {
                await window.mantis.bags.create('TW: ' + ph.name, ids);
                made += 1;
            } catch (e) { /* bag may already exist */ }
        }
        toast(made ? 'Wrote ' + made + ' commitment-phase bags.' : 'Could not write bags.');
    }

    async function drawFlow() {
        // Downstream edges aggregated into a star from each stage's most-committed
        // cell toward the cells it flows into, which is what the SDK draws.
        const byId = new Map(CELLS.map((c) => [c.twId, c]));
        let focus = null, mostQ = -1;
        for (const c of CELLS) if (c.downstream.length && c.q > mostQ) { mostQ = c.q; focus = c; }
        if (!focus) { toast('No downstream flow on this map.'); return; }
        const to = focus.downstream.map((d) => byId.get(d)).filter(Boolean).map((c) => c.id);
        if (!to.length) { toast('No downstream targets found.'); return; }
        try {
            await window.mantis.edges.set({ from: focus.id, to, predicate: 'flows to' });
            toast('Drew ' + to.length + ' downstream edges from the leading cell.');
        } catch (e) { toast('Could not draw flow: ' + e.message); }
    }

    function showCell(cell) {
        const el = $('twCell');
        if (!el) return;
        if (!cell) {
            el.className = 'cell';
            el.innerHTML = '<span class="hint-dot"></span>Click a cell on the map to place it on the committor and trace its path to commitment.';
            return;
        }
        // Trace a greedy path to commitment: follow the most-committed downstream
        // neighbour until q >= 0.9 or the chain ends.
        const byId = new Map(CELLS.map((c) => [c.twId, c]));
        const path = [cell];
        const seen = new Set([cell.twId]);
        let cur = cell;
        for (let i = 0; i < 40 && cur.q < 0.9; i += 1) {
            let next = null;
            for (const d of cur.downstream) {
                const c = byId.get(d);
                if (c && !seen.has(c.twId) && (!next || c.q > next.q)) next = c;
            }
            if (!next) break;
            path.push(next); seen.add(next.twId); cur = next;
        }
        el.className = 'cell filled';
        el.innerHTML = '<div class="cell-top">'
            + '<span class="cell-stage">' + escapeHtml(cell.stage) + '</span>'
            + '<span class="cell-phase" style="color:' + phaseColor(cell.phase) + '">' + escapeHtml(cell.phase) + '</span>'
            + '</div>'
            + '<div class="cell-grid">'
            + '<div><span class="cl">committor</span><span class="cv">' + fmt(cell.q, 3)
            + (cell.sd !== null ? ' <span class="pm">&plusmn;' + fmt(cell.sd, 2) + '</span>' : '') + '</span></div>'
            + '<div><span class="cl">free energy</span><span class="cv">' + (cell.energy === null ? '--' : fmt(cell.energy, 2) + ' ' + KT) + '</span></div>'
            + '<div><span class="cl">irreversibility</span><span class="cv">' + (cell.excess === null ? '--' : (cell.excess > 0 ? '+' : '') + fmt(cell.excess, 3)) + '</span></div>'
            + '<div><span class="cl">downstream</span><span class="cv">' + cell.downstream.length + '</span></div>'
            + '</div>'
            + '<p class="caption mono">' + path.map((c) => c.q.toFixed(2)).join(' &rarr; ') + '</p>';
        // Draw the traced path on the map.
        if (path.length > 1) {
            window.mantis.edges.set({ from: cell.id, to: path.slice(1).map((c) => c.id), predicate: 'commits via' })
                .catch(() => undefined);
        }
    }

    function wireInteractions() {
        $('twBags').addEventListener('click', makePhaseBags);
        $('twFlow').addEventListener('click', drawFlow);
        $('twPhases').addEventListener('click', (ev) => {
            const row = ev.target.closest('.phase');
            if (!row) return;
            const name = row.getAttribute('data-phase');
            const ids = CELLS.filter((c) => c.phase === name).map((c) => c.id);
            window.mantis.selection.set(ids).catch(() => undefined);
            window.mantis.bags.create('TW: ' + name, ids).then(() => toast('Selected ' + ids.length + ' ' + name + ' cells.')).catch(() => undefined);
        });
        // Redraw the field when the panel is resized.
        let rt = 0;
        window.addEventListener('resize', () => {
            clearTimeout(rt);
            rt = setTimeout(() => { if (CELLS.length) drawField($('twField'), CELLS); }, 150);
        });
    }

    async function onSelection() {
        try {
            const sel = await window.mantis.selection.get();
            const id = sel.currentSelection || (sel.selected && sel.selected[0]);
            const cell = CELLS.find((c) => c.id === id);
            if (cell) showCell(cell);
        } catch (e) { /* ignore */ }
    }

    function showEmpty(msg) {
        $('twBody').innerHTML = '<section class="card empty-card"><h2>No Thermodynamic Waddington export on this map</h2>'
            + '<p class="caption">' + escapeHtml(msg || 'This panel reads maps whose points carry tw_committor and the other tw_* fields. Open a map built from a Thermodynamic Waddington export.') + '</p></section>';
        $('twCount').textContent = '0 cells';
    }

    async function main() {
        const root = document.getElementById('root') || document.body;
        root.innerHTML = TEMPLATE;
        $('twCodewave').textContent = makeCodewave(12, 150);
        mountShader($('twShader'));

        if (!window.mantis || !window.mantis.maps) { showEmpty('The Mantis panel API is not available.'); return; }
        try {
            const active = await window.mantis.maps.getActive();
            MAP_ID = active && active.mapId;
            const points = await window.mantis.maps.getPoints(MAP_ID);
            CELLS = (points || []).map(cellOf).filter(Boolean);
        } catch (e) {
            showEmpty('Could not read the active map: ' + e.message);
            return;
        }
        if (!CELLS.length) { showEmpty(); return; }
        renderSummary(CELLS);
        wireInteractions();

        try {
            await window.mantis.events.subscribe('selection.changed', {}, onSelection);
            onSelection();
        } catch (e) { /* selection events optional */ }
        window.__TW_CELLS__ = CELLS.length;
    }

    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', () => { main().catch((e) => showEmpty(e.message)); });
    } else {
        main().catch((e) => showEmpty(e.message));
    }
})();
