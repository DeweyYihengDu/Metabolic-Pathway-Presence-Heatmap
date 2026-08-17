"""Neighbour calibration at scale: does a closer training genome calibrate better?"""
import csv
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "benchmarks/threshold_audit")
sys.path.insert(0, "benchmarks/annotation")
from calibration_transfer import brier, ece, fit_logistic, predict
from compare_annotation import load_ground_truth, load_protein_id_map

root = Path("panel/data_gb")
_tax_path = (root / "taxonomy_panel.tsv" if (root / "taxonomy_panel.tsv").exists()
             else Path("panel/calibration_panel.tsv"))
with open(_tax_path, encoding="utf-8") as _fh:
    tax = {r["org"]: r for r in csv.DictReader(_fh, delimiter="\t")}

data = {}
for org in sorted(tax):
    hits = root/"subthreshold"/f"subthreshold_{org}.tsv"
    gt = root/"ground_truth"/f"link_ko_{org}.tsv"
    idm = root/"ground_truth"/f"protein_id_map_{org}.tsv"
    if not (hits.exists() and gt.exists() and idm.exists()):
        continue
    df = pd.read_csv(hits, sep="\t")
    df["gene"] = df["gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    m = load_protein_id_map(idm, org)
    truth = load_ground_truth(gt, org, m)
    df = df[df["gene"].isin(set(m.values()))]
    tp = {(g,k) for g,kos in truth.items() for k in kos}
    acc = df[df["delta"] >= 0]
    if len(acc) < 200:
        continue
    y = np.array([(g,k) in tp for g,k in zip(acc["gene"], acc["ko"])], float)
    data[org] = (acc["delta"].to_numpy(float), y)
print(f"{len(data)} organisms usable", flush=True)

def tier(a,b):
    ta, tb = tax[a], tax[b]
    if ta["genus"] and ta["genus"]==tb["genus"]: return "1_same_genus"
    if ta["clade"]==tb["clade"]: return "2_same_clade"
    return "3_same_domain"

rows=[]
orgs=sorted(data)
for train, test in itertools.permutations(orgs, 2):
    d_tr,y_tr = data[train]; d_te,y_te = data[test]
    a,b = fit_logistic(d_tr,y_tr)
    pm = predict(d_te,a,b); pc = np.full_like(y_te, float(y_tr.mean()))
    rows.append({"train":train,"test":test,"tier":tier(train,test),
                 "ece_model":ece(pm,y_te),"ece_constant":ece(pc,y_te),
                 "brier_model":brier(pm,y_te),"brier_constant":brier(pc,y_te)})
df=pd.DataFrame(rows)
df["ece_improvement"]=df.ece_constant-df.ece_model
df["brier_improvement"]=df.brier_constant-df.brier_model
df.to_csv("panel/neighbour_scale.csv", index=False)
print("\n=== mean by relatedness tier ===")
print(df.groupby("tier")[["ece_model","ece_constant","brier_model","brier_improvement"]]
        .agg(["mean","count"]).round(4).to_string())
