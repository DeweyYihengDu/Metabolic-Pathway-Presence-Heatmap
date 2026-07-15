"""Build a self-contained interactive HTML report from an MPPH output folder.

The report embeds the matrix, feature table, QC table and manifest as JSON and
renders a colour-coded, searchable heatmap with hover tooltips -- no external
assets, so it opens offline from a single file.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

_TEMPLATE = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MPPH report · {slug}</title>
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
<header><h1>MPPH report &middot; {slug}</h1>
<div class="sub">{subtitle}</div></header>
<main>
<div class="tabs">
<button data-tab="heatmap" class="active">Heatmap</button>
<button data-tab="qc">QC ({n_org} organisms)</button>
<button data-tab="manifest">Manifest</button>
</div>
<section id="heatmap">
<div class="controls">
<input type="search" id="filter" placeholder="Filter features by id or name…">
<span class="legend" id="legend"></span>
</div>
<div class="scroll"><table id="grid"></table></div>
</section>
<section id="qc" class="hidden"><div class="scroll"><table class="qc" id="qctab"></table></div></section>
<section id="manifest" class="hidden"><div class="scroll"><table class="mani" id="manitab"></table></div></section>
</main>
<script>
const DATA = {data_json};
const mode = DATA.mode;
function color(v){{
  if(v<=0) return '#eceef1';
  if(mode==='completeness'){{
    const t=Math.min(1,v); const l=Math.round(240-150*t), g=Math.round(238-90*t), b=Math.round(241-70*t);
    return 'rgb('+Math.round(238-200*t)+','+g+','+Math.round(251-160*t)+')';
  }}
  return '#2a78d6';
}}
function buildGrid(){{
  const t=document.getElementById('grid');
  const cols=DATA.features, rows=DATA.organisms, M=DATA.matrix;
  let h='<thead><tr><th class="rowh">organism</th>';
  cols.forEach((c,j)=>{{h+='<th data-col="'+j+'" title="'+c.id+' '+c.name+'">'+c.id+'</th>';}});
  h+='</tr></thead><tbody>';
  rows.forEach((r,i)=>{{
    h+='<tr><td class="rowh">'+r+'</td>';
    cols.forEach((c,j)=>{{const v=M[i][j];
      h+='<td class="cell" data-col="'+j+'" style="background:'+color(v)+'" title="'+r+' × '+c.id+' '+c.name+' = '+v.toFixed(2)+'"></td>';}});
    h+='</tr>';
  }});
  t.innerHTML=h+'</tbody>';
  document.getElementById('legend').textContent =
    mode==='completeness' ? 'cell colour = module completeness (0→1)' : 'blue = present, grey = absent';
}}
function buildTable(id,obj){{
  const rows=Object.entries(obj).map(([k,v])=>'<tr><th>'+k+'</th><td>'+
    (typeof v==='object'?JSON.stringify(v):v)+'</td></tr>').join('');
  document.getElementById(id).innerHTML=rows;
}}
function buildQC(){{
  const q=DATA.qc; if(!q||!q.length){{return;}}
  const keys=Object.keys(q[0]);
  let h='<thead><tr>'+keys.map(k=>'<th>'+k+'</th>').join('')+'</tr></thead><tbody>';
  q.forEach(r=>{{h+='<tr>'+keys.map(k=>'<td>'+r[k]+'</td>').join('')+'</tr>';}});
  document.getElementById('qctab').innerHTML=h+'</tbody>';
}}
buildGrid(); buildQC(); buildTable('manitab',DATA.manifest);
document.querySelectorAll('.tabs button').forEach(b=>b.onclick=()=>{{
  document.querySelectorAll('.tabs button').forEach(x=>x.classList.remove('active'));
  b.classList.add('active');
  ['heatmap','qc','manifest'].forEach(s=>document.getElementById(s).classList.add('hidden'));
  document.getElementById(b.dataset.tab).classList.remove('hidden');
}});
document.getElementById('filter').oninput=e=>{{
  const q=e.target.value.toLowerCase();
  DATA.features.forEach((c,j)=>{{
    const show=!q||c.id.toLowerCase().includes(q)||c.name.toLowerCase().includes(q);
    document.querySelectorAll('[data-col="'+j+'"]').forEach(el=>el.style.display=show?'':'none');
  }});
}};
</script></body></html>"""


def build_report(outdir: str | Path, slug: str) -> Path:
    """Read the ``<slug>_*`` outputs in ``outdir`` and write ``<slug>_report.html``."""
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
    data = {
        "mode": manifest.get("mode", "presence"),
        "organisms": [str(i) for i in matrix.index],
        "features": feats,
        "matrix": matrix.round(4).to_numpy().tolist(),
        "qc": qc,
        "manifest": manifest,
    }
    subtitle = (f"{manifest.get('mode', '')} &middot; "
                f"{matrix.shape[0]} organisms &times; {matrix.shape[1]} features "
                f"&middot; MPPH {manifest.get('mpph_version', '')}")
    html = _TEMPLATE.format(slug=slug, subtitle=subtitle, n_org=matrix.shape[0],
                            data_json=json.dumps(data))
    out = outdir / f"{slug}_report.html"
    out.write_text(html, encoding="utf-8")
    return out
