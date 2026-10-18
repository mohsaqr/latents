"""Read-only CSV/output checks; writes a log only when stdout is redirected."""
import csv
import itertools
from pathlib import Path
import math
import random
import re

folder = Path("equivalence/starts-study")
rows = list(csv.DictReader((folder / "blocks.csv").open()))
comparison = list(csv.DictReader((folder / "mplus-comparison.csv").open()))
cell_rows = list(csv.DictReader((folder / "codex/cell-diagnostics.csv").open()))
cells = {(r["dataset"], int(r["k"])): r for r in cell_rows}
outputs = list((folder / "mplus").glob("*.out"))
finished = 0
for path in outputs:
    text = path.read_text()
    if "THE MODEL ESTIMATION TERMINATED NORMALLY" not in text:
        continue
    finished += 1
    ranked = text.split("Final stage loglikelihood values at local maxima", 1)[1]
    values = []
    for line in ranked.splitlines()[2:]:
        if not line.strip():
            break
        values.append(float(line.split()[0]))
    h0 = float(re.search(r"^\s+H0 Value\s+(-?[\d.]+)", text, re.M)[1])
    assert values and values == sorted(values, reverse=True) and values[0] == h0
print("OUTPUTS", len(outputs), "NORMAL", finished, "UNFINISHED", len(outputs) - finished)
for r in comparison:
    path = folder / "mplus" / (r["dataset"] + "-k" + r["k"] + ".out")
    text = path.read_text()
    h0 = float(re.search(r"^\s+H0 Value\s+(-?[\d.]+)", text, re.M)[1])
    assert h0 == float(r["mplus_best"])
    assert abs(float(cells[r["dataset"], int(r["k"])]["reference"]) - float(r["reference"])) < 1e-8
    assert abs(float(r["reference"]) - h0 - float(r["latents_minus_mplus"])) < 1e-8
print("COMPARISON ROWS VERIFIED", len(comparison))
lpa_over = [r for r in comparison if r["family"] == "lpa" and r["relative_k"] == "K+1"]
zero = nonpd = 0
for r in lpa_over:
    text = (folder / "mplus" / (r["dataset"] + "-k" + r["k"] + ".out")).read_text()
    model = text.split("MODEL RESULTS", 1)[1].split("QUALITY OF NUMERICAL RESULTS", 1)[0]
    variance_values = []
    reading = False
    for line in model.splitlines():
        if line.strip() == "Variances":
            reading = True
            continue
        if reading:
            match = re.match(r"\s+Y[1-4]\s+(-?[\d.]+)", line)
            if match:
                variance_values.append(float(match[1]))
            elif line.strip():
                reading = False
    assert len(variance_values) == 4 * int(r["k"])
    zero += min(variance_values) == 0
    nonpd += "NON-POSITIVE DEFINITE" in text
print("LPA K+1 MATCHED", len(lpa_over), "DISPLAYED ZERO VARIANCE", zero, "NONPD", nonpd)

# Shuffle each k's block labels independently, retaining every likelihood.
rng = random.Random(20261018)
for budget in (10, 20):
    selected = [r for r in rows if int(r["n_starts"]) == budget]
    groups = {}
    for r in selected:
        groups.setdefault(r["dataset"], {}).setdefault(int(r["k"]), []).append(r)
    records = []
    for _ in range(1000):
        fewer = more = 0
        for dataset, ks in groups.items():
            ordered = sorted(ks)
            bics = []
            ref = []
            for k in ordered:
                z = ks[k]
                p = int(z[0]["n_parameters"])
                penalty = p * math.log(int(z[0]["n"]))
                ref.append(-2 * float(z[0]["reference"]) + penalty)
                values = [-2 * float(r["log_likelihood"]) + penalty for r in z]
                rng.shuffle(values)
                bics.append(values)
            reference = ordered[min(range(len(ordered)), key=ref.__getitem__)]
            for run in zip(*bics):
                chosen = ordered[min(range(len(ordered)), key=run.__getitem__)]
                fewer += chosen < reference
                more += chosen > reference
        records.append((fewer, more))
    print("PAIRING", budget, "1000 SHUFFLES FEWER RANGE", min(x[0] for x in records), max(x[0] for x in records),
          "MORE RANGE", min(x[1] for x in records), max(x[1] for x in records))
print("ALL CHECKS PASSED")
