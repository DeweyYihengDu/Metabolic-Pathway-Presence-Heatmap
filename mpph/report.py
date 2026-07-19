"""Build a self-contained interactive HTML report from an MPPH output folder.

Security: organism/feature/metadata strings are user-controlled, so the payload
is embedded in a ``<script type="application/json">`` block with ``<`` escaped
(no ``</script>`` breakout) and every cell is written with ``textContent`` via
DOM construction -- never ``innerHTML`` on user data. A CSP is set as defence in
depth.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd

_TEMPLATE = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; \
script-src 'unsafe-inline'; style-src 'unsafe-inline'; img-src data:;">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MPPH report</title>
<style>
:root{{--ink:#0b0b0b;--muted:#898781;--line:#e1e0d9;--surface:#fcfcfb}}
*{{box-sizing:border-box}}
body{{font-family:system-ui,Segoe UI,Arial,sans-serif;margin:0;color:var(--ink);
background:var(--surface)}}
header{{padding:18px 22px;border-bottom:1px solid var(--line)}}
h1{{font-size:19px;margin:0 0 4px}} .sub{{color:var(--muted);font-size:13px}}
main{{padding:16px 22px}}
.controls{{margin:10px 0;display:flex;gap:10px;flex-wrap:wrap;align-items:center}}
input[type=search]{{padding:6px 10px;border:1px solid var(--line);border-radius:6px;
font-size:13px;min-width:220px}}
.tabs button{{border:1px solid var(--line);background:#fff;padding:6px 12px;
border-radius:6px;cursor:pointer;font-size:13px}}
.tabs button.active{{background:#2a78d6;color:#fff;border-color:#2a78d6}}
.scroll{{overflow:auto;max-height:70vh;border:1px solid var(--line);border-radius:8px}}
table{{border-collapse:collapse;font-size:12px}}
th,td{{border:1px solid #f0efec;padding:0}}
th.rowh,td.rowh{{padding:3px 8px;position:sticky;left:0;background:#fff;
white-space:nowrap;font-style:italic;text-align:left;z-index:1}}
thead th{{position:sticky;top:0;background:#fff;writing-mode:vertical-rl;
transform:rotate(180deg);padding:4px 2px;font-weight:500;color:var(--muted)}}
td.cell{{width:16px;height:16px}}
.qc td,.qc th,.mani td,.mani th{{border:1px solid var(--line);padding:4px 8px;
font-size:12px;text-align:left}}
.hidden{{display:none}}
.legend{{font-size:12px;color:var(--muted);margin-left:auto}}
</style></head><body>
<header><h1 id="title"></h1><div class="sub" id="subtitle"></div></header>
<main>
<div class="tabs">
<button data-tab="heatmap" class="active">Heatmap</button>
<button data-tab="qc">QC</button>
<button data-tab="manifest">Manifest</button>
</div>
<section id="heatmap">
<div class="controls">
<input type="search" id="filter" placeholder="Filter features by id or name...">
<span class="legend" id="legend"></span>
</div>
<div class="scroll"><table id="grid"></table></div>
</section>
<section id="qc" class="hidden"><div class="scroll"><table class="qc" id="qctab"></table></div></section>
<section id="manifest" class="hidden"><div class="scroll"><table class="mani" id="manitab"></table></div></section>
</main>
<script id="mpph-data" type="application/json">{data_json}</script>
<script>
const DATA = JSON.parse(document.getElementById('mpph-data').textContent);
const mode = DATA.mode;
function el(tag, text){{
  const e = document.createElement(tag);
  if (text !== undefined && text !== null) e.textContent = String(text);
  return e;
}}
function color(v){{
  if (v === null || !(v > 0)) return '#eceef1';
  if (mode === 'completeness'){{
    const t = Math.min(1, v);
    return 'rgb(' + Math.round(238-200*t) + ',' + Math.round(238-90*t) + ','
      + Math.round(251-160*t) + ')';
  }}
  return '#2a78d6';
}}
function buildGrid(){{
  const table = document.getElementById('grid');
  if (DATA.grid_omitted){{
    document.getElementById('filter').style.display = 'none';
    document.getElementById('legend').textContent =
      DATA.organisms.length + ' organisms x ' + DATA.features.length +
      ' features (' + (DATA.organisms.length * DATA.features.length).toLocaleString() +
      ' cells) is too large for this interactive view (over ' +
      DATA.max_cells.toLocaleString() + ' cells) -- a browser table that size ' +
      'can hang the page. Full data: ' + DATA.slug + '_matrix.csv / ' +
      DATA.slug + '_ordered_matrix.csv, and the heatmap figure(s).';
    return;
  }}
  const cols = DATA.features, rows = DATA.organisms, M = DATA.matrix;
  const thead = el('thead'), hr = el('tr');
  hr.appendChild(Object.assign(el('th','organism'), {{className:'rowh'}}));
  cols.forEach((c, j) => {{
    const th = el('th', c.id);
    th.dataset.col = j; th.title = c.id + ' ' + c.name;
    hr.appendChild(th);
  }});
  thead.appendChild(hr); table.appendChild(thead);
  const tb = el('tbody');
  rows.forEach((r, i) => {{
    const tr = el('tr');
    tr.appendChild(Object.assign(el('td', r), {{className:'rowh'}}));
    cols.forEach((c, j) => {{
      const v = M[i][j];
      const td = el('td');
      td.className = 'cell'; td.dataset.col = j;
      td.style.background = color(v);
      td.title = r + ' x ' + c.id + ' ' + c.name + ' = '
        + (v === null ? 'unknown' : v.toFixed(2));
      tr.appendChild(td);
    }});
    tb.appendChild(tr);
  }});
  table.appendChild(tb);
  document.getElementById('legend').textContent = mode === 'completeness'
    ? 'cell colour = module completeness (0 to 1)'
    : 'blue = present, grey = absent';
}}
function buildRecords(id, records){{
  const t = document.getElementById(id);
  if (!records || !records.length) return;
  const keys = Object.keys(records[0]);
  const thead = el('thead'), hr = el('tr');
  keys.forEach(k => hr.appendChild(el('th', k)));
  thead.appendChild(hr); t.appendChild(thead);
  const tb = el('tbody');
  records.forEach(rec => {{
    const tr = el('tr');
    keys.forEach(k => tr.appendChild(el('td', rec[k])));
    tb.appendChild(tr);
  }});
  t.appendChild(tb);
}}
function buildPairs(id, obj){{
  const t = document.getElementById(id), tb = el('tbody');
  Object.entries(obj).forEach(([k, v]) => {{
    const tr = el('tr');
    tr.appendChild(el('th', k));
    tr.appendChild(el('td', typeof v === 'object' ? JSON.stringify(v) : v));
    tb.appendChild(tr);
  }});
  t.appendChild(tb);
}}
document.getElementById('title').textContent = 'MPPH report - ' + DATA.slug;
document.getElementById('subtitle').textContent = DATA.subtitle;
buildGrid(); buildRecords('qctab', DATA.qc); buildPairs('manitab', DATA.manifest);
document.querySelectorAll('.tabs button').forEach(b => b.onclick = () => {{
  document.querySelectorAll('.tabs button').forEach(x => x.classList.remove('active'));
  b.classList.add('active');
  ['heatmap','qc','manifest'].forEach(s =>
    document.getElementById(s).classList.add('hidden'));
  document.getElementById(b.dataset.tab).classList.remove('hidden');
}});
document.getElementById('filter').oninput = e => {{
  const q = e.target.value.toLowerCase();
  DATA.features.forEach((c, j) => {{
    const show = !q || c.id.toLowerCase().includes(q)
      || c.name.toLowerCase().includes(q);
    document.querySelectorAll('[data-col="' + j + '"]').forEach(
      x => x.style.display = show ? '' : 'none');
  }});
}};
</script></body></html>"""


def _safe_json(payload: dict) -> str:
    """JSON for a <script type="application/json"> block.

    ``</script>`` (and any other tag) cannot break out because ``<`` is escaped;
    NaN/Infinity are rejected by allow_nan and converted to null beforehand.
    """
    return (json.dumps(payload, ensure_ascii=False, allow_nan=False)
            .replace("<", "\\u003c")
            .replace(" ", "\\u2028")
            .replace(" ", "\\u2029"))


def _clean(value):
    """Convert NaN/Infinity to None so the payload is valid JSON."""
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def build_report(outdir: str | Path, slug: str, *, max_cells: int = 50_000) -> Path:
    """Read the ``<slug>_*`` outputs in ``outdir`` and write ``<slug>_report.html``.

    A matrix with more than ``max_cells`` organisms x features is not embedded
    as an interactive table -- an ``organisms x features`` browser table that
    large (each cell its own styled DOM node) can hang the page. The QC and
    manifest tabs, and a summary in place of the grid, are still built; the
    full data is always in ``<slug>_matrix.csv`` / ``_ordered_matrix.csv`` and
    the heatmap figure(s) regardless.
    """
    outdir = Path(outdir)
    matrix = pd.read_csv(outdir / f"{slug}_matrix.csv", index_col=0)
    matrix.columns = matrix.columns.astype(str)
    features = pd.read_csv(outdir / f"{slug}_features.csv", dtype={"feature_id": str})
    manifest = json.loads((outdir / f"{slug}_manifest.json").read_text("utf-8"))
    qc_path = outdir / f"{slug}_qc.csv"
    qc = pd.read_csv(qc_path).to_dict("records") if qc_path.exists() else []

    name_by_id = dict(zip(features["feature_id"].astype(str), features["name"]))
    feats = [{"id": str(c), "name": str(name_by_id.get(str(c), c))}
             for c in matrix.columns]
    n_cells = matrix.shape[0] * matrix.shape[1]
    grid_omitted = n_cells > max_cells
    values = [] if grid_omitted else [
        [_clean(v) for v in row] for row in matrix.round(4).to_numpy().tolist()]
    subtitle = (f"{manifest.get('mode', '')} - {matrix.shape[0]} organisms x "
                f"{matrix.shape[1]} features - MPPH "
                f"{manifest.get('mpph_version', '')}")
    data = {
        "slug": slug,
        "subtitle": subtitle,
        "mode": manifest.get("mode", "presence"),
        "grid_omitted": grid_omitted,
        "max_cells": max_cells,
        "organisms": [str(i) for i in matrix.index],
        "features": feats,
        "matrix": values,
        "qc": [{k: _clean(v) for k, v in rec.items()} for rec in qc],
        "manifest": manifest,
    }
    html = _TEMPLATE.format(data_json=_safe_json(data))
    out = outdir / f"{slug}_report.html"
    out.write_text(html, encoding="utf-8")
    return out
