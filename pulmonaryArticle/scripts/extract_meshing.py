#!/usr/bin/env python3
"""Collect the meshing-robustness data for Section 3.1 into one CSV.

Reads a directory of fsiScript meshing runs (one meshing.sh run per geometry,
as made by the cohort test) and the matching surface probe results, and
writes one row per geometry to data/meshing_ext8.csv:

  * outcome: clean (solid passes the standard OpenFOAM checks), orientation
    (only incorrectly oriented faces / quality failures, no negative volumes),
    negative (negative-volume cells), or failed (no volume meshes);
  * failure cause for failed runs, classified from the meshing log;
  * fluid and solid cell counts and mesh-quality measures from checkMesh;
  * the share of surface points moved by the windowed-sinc smoothing and by
    the saddle rounding, and the change in lumen volume (surface_probe.py);
  * the share of wall points whose extrusion normal was smoothed, from
    varExtrudeMesh's log.

Usage:
  extract_meshing.py RUN_DIR PROBE_DIR [OUTPUT_CSV]

RUN_DIR holds results/<geometry>.json (written by the cohort test's parse
step), cases/log.<geometry> and cases/caseFiles/<geometry>/; PROBE_DIR holds
<geometry>.json from surface_probe.py.
"""
import csv
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
run_dir = Path(sys.argv[1])
probe_dir = Path(sys.argv[2])
output = Path(sys.argv[3]) if len(sys.argv) > 3 else HERE.parent / "data" / "meshing_ext8.csv"

FIELDS = [
    "geometry", "species", "outcome", "failure_cause", "openings",
    "fluid_cells", "solid_cells",
    "fluid_max_nonortho", "fluid_avg_nonortho", "fluid_max_skewness", "fluid_max_aspect",
    "solid_max_nonortho", "solid_avg_nonortho", "solid_max_skewness", "solid_max_aspect",
    "solid_negative_cells", "solid_wrong_faces",
    "taubin_moved_fraction", "taubin_mean_displacement_mm",
    "saddle_moved_fraction", "saddle_mean_displacement_moved_mm", "saddle_max_displacement_mm",
    "normal_smoothed_fraction", "normal_rounds", "normal_faces_folding_before",
    "volume_change_taubin", "volume_change_saddle", "volume_change_total",
    "runtime_s",
]


def classify_failure(result, log_text):
    if result.get("stage") == "outlet shut" or "meshed shut" in log_text:
        return "outlet closed"
    error = result.get("error") or ""
    m = re.search(r"\((\d+) edge\(s\) on the wall, (\d+) on the caps\)", log_text)
    if m:
        wall, caps = int(m.group(1)), int(m.group(2))
        return "touching openings" if caps > 0 and wall == 0 else "non-manifold wall in scan"
    if "not a simple closed loop" in log_text:
        return "overlapping duplicated patch in scan" if "duplicate triangle" in log_text \
            else "rim not a simple loop"
    if "step failed: Meshing the fluid" in log_text:
        return "fluid mesher failed"
    return "other: " + error[:80]


def adaptive_stats(case):
    log = case / "varExtrudeFunc" / "log.varExtrudeMesh"
    if not log.exists():
        return {}
    text = log.read_text(errors="replace")
    out = {}
    m = re.search(r"points whose normal was smoothed: (\d+) of (\d+)", text)
    if m:
        out["normal_smoothed_fraction"] = int(m.group(1)) / int(m.group(2))
    m = re.search(r"rounds: (\d+)", text)
    if m:
        out["normal_rounds"] = int(m.group(1))
    m = re.search(r"faces that would turn over before: (\d+)", text)
    if m:
        out["normal_faces_folding_before"] = int(m.group(1))
    return out


rows = []
for path in sorted((run_dir / "results").glob("*.json")):
    r = json.loads(path.read_text())
    g = r["geometry"]
    species = "porcine" if g.startswith("p") else "human"
    log_file = run_dir / "cases" / f"log.{g}"
    log_text = log_file.read_text(errors="replace") if log_file.exists() else ""
    row = {k: "" for k in FIELDS}
    row.update(geometry=g, species=species, openings=r.get("profiles") or "",
               runtime_s=r.get("seconds", ""))
    if r["rc"] != 0:
        row["outcome"] = "failed"
        row["failure_cause"] = classify_failure(r, log_text)
    else:
        neg = r.get("solid_negative") or 0
        row["outcome"] = "clean" if r.get("solid_status") == "ok" else ("negative" if neg else "orientation")
        for k in ("fluid_cells", "solid_cells", "fluid_max_nonortho", "fluid_avg_nonortho",
                  "fluid_max_skewness", "fluid_max_aspect", "solid_max_nonortho",
                  "solid_avg_nonortho", "solid_max_skewness", "solid_max_aspect"):
            row[k] = r.get(k, "")
        row["solid_negative_cells"] = neg
        row["solid_wrong_faces"] = r.get("solid_wrong_faces") or 0
        row.update(adaptive_stats(run_dir / "cases" / "caseFiles" / g))
    probe = probe_dir / f"{g}.json"
    if probe.exists():
        p = json.loads(probe.read_text())
        row.update(
            taubin_moved_fraction=p["taubin"]["fraction_moved"],
            taubin_mean_displacement_mm=p["taubin"]["mean_displacement_all"],
            saddle_moved_fraction=p["saddle"]["fraction_moved"],
            saddle_mean_displacement_moved_mm=p["saddle"]["mean_displacement_moved"],
            saddle_max_displacement_mm=p["saddle"]["max_displacement"],
            volume_change_taubin=p["volume"]["change_taubin"],
            volume_change_saddle=p["volume"]["change_saddle"],
            volume_change_total=p["volume"]["change_total"],
        )
    rows.append(row)

output.parent.mkdir(parents=True, exist_ok=True)
with output.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, FIELDS)
    writer.writeheader()
    writer.writerows(rows)
print(f"Wrote {len(rows)} geometries to {output}")
