"""Recover compact real Airway figure from the previously qualified artifact."""
import csv, hashlib, io, json, math, sys, zipfile
from pathlib import Path
if len(sys.argv)!=3: raise SystemExit("usage: python scripts/extract_airway_volcano.py archive.zip output.json")
p=Path(sys.argv[1])
assert hashlib.sha256(p.read_bytes()).hexdigest()=="c9fe3e1ae37892b0c263e7702171225e64f0659f7ef7ca4df47d645b888bcc1d"
with zipfile.ZipFile(p) as z:
    report=json.loads(z.read("paired_airway_original_response.json"))["methods"]["transcriptomics_differential_deseq2"]
assert report["referenceTableComplete"] is True
rows=list(csv.DictReader(io.StringIO(report["referenceTableCsv"])))
assert len(rows)==17199 and len({r["feature"] for r in rows})==17199
bins=[[0]*16 for _ in range(12)]
valid=[]
for row in rows:
    x,q=float(row["log2FoldChange"]),float(row["padj"])
    assert math.isfinite(x) and -6<=x<10 and math.isfinite(q) and 0<q<=1
    valid.append((x,q))
    bins[min(11,math.floor(-math.log10(q)))][math.floor(x)+6]+=1
q=[.01,.05,.1]
fx=[0,.5,1,2]
counts=[[sum(a<=t and abs(x)>=f for x,a in valid) for f in fx] for t in q]
assert sum(map(sum,bins))==17199 and [c[0] for c in counts]==[2940,4100,4947]
marker=next(row for row in rows if row["feature"]=="ENSG00000103196")
data=dict(n=len(valid),bins=bins,q=q,fx=fx,counts=counts,marker=dict(feature="CRISPLD2",log2FC=float(marker["log2FoldChange"]),q=float(marker["padj"])))
Path(sys.argv[2]).write_text(json.dumps(data,ensure_ascii=False,separators=(",",":"))+"\n",encoding="utf-8")
