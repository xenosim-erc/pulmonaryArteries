#!/usr/bin/env python3
"""Section 3.1 (pipeline robustness): table and numbers from the meshing data.

Reads data/meshing_ext8.csv (made by extract_meshing.py) and writes
  tables/meshing.tex         the summary table, \\input in results.tex
  tables/meshingNumbers.tex  \\newcommand macros for every number the text
                             quotes, \\input in results.tex

Statistics of mesh size, quality and surface modification are over the
geometries whose solid mesh passed the standard checks ("clean").
"""
import csv
import statistics as st
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
rows = list(csv.DictReader((ROOT / "data" / "meshing_ext8.csv").open()))
GROUPS = {"Human": [r for r in rows if r["species"] == "human"],
          "Porcine": [r for r in rows if r["species"] == "porcine"],
          "All": rows}
CAUSES = [("outlet closed", "OutletClosed", "small outlet closed"),
          ("touching openings", "TouchingOpenings", "openings touching (scan)"),
          ("non-manifold wall in scan", "NonManifoldWall", "non-manifold wall (scan)"),
          ("overlapping duplicated patch in scan", "DuplicatedPatch", "duplicated wall patch (scan)")]


def values(group, key, outcome="clean"):
    return [float(r[key]) for r in group if r["outcome"] == outcome and r[key] != ""]


def quartiles(v):
    q = st.quantiles(v, n=4, method="inclusive")
    return q[0], st.median(v), q[2]


macros = {}


def macro(name, value):
    macros[name] = value


table = []
fmt_int = lambda x: f"{x:d}"


def row(label, cells, indent=False):
    table.append(("\\quad " if indent else "") + label + " & " + " & ".join(cells) + r" \\")


def count_row(label, key_fn, name, indent=False):
    cells = []
    for g, group in GROUPS.items():
        n = sum(1 for r in group if key_fn(r))
        cells.append(fmt_int(n))
        macro(f"Mesh{name}{g}", fmt_int(n))
        macro(f"Mesh{name}Pct{g}", f"{100 * n / len(group):.0f}")
    row(label, cells, indent)


def stat_row(label, key, name, scale=1.0, digits=0, kind="iqr"):
    cells = []
    for g, group in GROUPS.items():
        v = [x * scale for x in values(group, key)]
        if not v:
            cells.append("--")
            continue
        q1, med, q3 = quartiles(v)
        lo, hi = min(v), max(v)
        f = lambda x: f"{x:.{digits}f}".replace("-", "$-$")
        if kind == "iqr":
            cells.append(f"{f(med)} [{f(q1)}, {f(q3)}]")
        elif kind == "max":
            cells.append(f"{f(med)} [{f(hi)}]")
        else:  # range
            cells.append(f"{f(med)} [{f(lo)}, {f(hi)}]")
        macro(f"{name}Median{g}", f(med))
        macro(f"{name}QOne{g}", f(q1))
        macro(f"{name}QThree{g}", f(q3))
        macro(f"{name}Min{g}", f(lo))
        macro(f"{name}Max{g}", f(hi))
    row(label, cells)


def section(title):
    table.append(r"\addlinespace" if table else "")
    table.append(r"\multicolumn{4}{l}{\textit{" + title + r"}} \\")


section("Outcome (geometries)")
count_row("Geometries", lambda r: True, "Geometries")
count_row("Volume meshes generated", lambda r: r["outcome"] != "failed", "Meshed")
count_row("Solid passes mesh checks", lambda r: r["outcome"] == "clean", "Clean", indent=True)
count_row("Oriented-face errors only", lambda r: r["outcome"] == "orientation", "Orientation", indent=True)
count_row("Negative-volume cells", lambda r: r["outcome"] == "negative", "Negative", indent=True)
count_row("No volume meshes", lambda r: r["outcome"] == "failed", "Failed")
for cause, name, label in CAUSES:
    count_row(label[0].upper() + label[1:], lambda r, c=cause: r["failure_cause"] == c, name, indent=True)
for g, group in GROUPS.items():
    n = sum(1 for r in group if r["failure_cause"] in [c for c, _, _ in CAUSES[1:]])
    macro(f"MeshScanDefects{g}", fmt_int(n))
other = [r for r in rows if r["outcome"] == "failed" and r["failure_cause"] not in [c for c, _, _ in CAUSES]]
if other:
    raise SystemExit(f"unclassified failures: {[(r['geometry'], r['failure_cause']) for r in other]}")

section(r"Mesh size ($\times 10^{3}$ cells)")
stat_row("Fluid cells", "fluid_cells", "FluidCells", scale=1e-3)
stat_row("Solid cells", "solid_cells", "SolidCells", scale=1e-3)

section("Mesh quality")
stat_row(r"Fluid max.\ non-orthogonality (\si{\degree})", "fluid_max_nonortho", "FluidNonOrtho", kind="max")
stat_row("Fluid max.\\ skewness", "fluid_max_skewness", "FluidSkew", digits=2, kind="max")
stat_row(r"Solid max.\ non-orthogonality (\si{\degree})", "solid_max_nonortho", "SolidNonOrtho", kind="max")
stat_row("Solid max.\\ skewness", "solid_max_skewness", "SolidSkew", digits=2, kind="max")

section(r"Surface, extrusion and lumen volume changes (\%)")
stat_row("Points moved: smoothing", "taubin_moved_fraction", "TaubinMoved", scale=100, digits=1, kind="range")
stat_row("Points moved: saddle rounding", "saddle_moved_fraction", "SaddleMoved", scale=100, digits=1, kind="range")
stat_row("Wall points: normal smoothed", "normal_smoothed_fraction", "NormalSmoothed", scale=100, digits=1, kind="range")
stat_row("Volume change: smoothing", "volume_change_taubin", "VolTaubin", scale=100, digits=2, kind="range")
stat_row("Volume change: saddle rounding", "volume_change_saddle", "VolSaddle", scale=100, digits=2, kind="range")
stat_row("Volume change: total", "volume_change_total", "VolTotal", scale=100, digits=2, kind="range")

# Further numbers the text quotes
clean = [r for r in rows if r["outcome"] == "clean"]
macro("TaubinDispMedianAll", f"{st.median(values(rows, 'taubin_mean_displacement_mm')):.2f}")
macro("SaddleDispMedianAll", f"{st.median(values(rows, 'saddle_mean_displacement_moved_mm')):.2f}")
macro("SaddleDispMaxAll", f"{max(values(rows, 'saddle_max_displacement_mm')):.1f}")
neg = [r for r in rows if r["outcome"] == "negative"]
macro("MeshNegativeCellsMax", fmt_int(max((int(float(r["solid_negative_cells"])) for r in neg), default=0)))
macro("MeshNegativeCellsTotal", fmt_int(sum(int(float(r["solid_negative_cells"])) for r in neg)))
macro("MeshRuntimeMedianMin", f"{st.median(float(r['runtime_s']) for r in rows) / 60:.0f}")
macro("MeshRuntimeMaxMin", f"{max(float(r['runtime_s']) for r in rows) / 60:.0f}")
macro("MeshOpeningsMedianHuman", f"{st.median(int(r['openings']) for r in GROUPS['Human'] if r['openings']):.0f}")
macro("MeshOpeningsMedianPorcine", f"{st.median(int(r['openings']) for r in GROUPS['Porcine'] if r['openings']):.0f}")

tables = ROOT / "tables"
tables.mkdir(exist_ok=True)
(tables / "meshing.tex").write_text("\n".join([
    "% Generated by scripts/meshing_results.py from data/meshing_ext8.csv; do not edit.",
    r"\begin{table}[tbp]",
    r"\centering",
    r"\footnotesize",
    r"\setlength{\tabcolsep}{3pt}",
    r"\caption{Outcome of the automated meshing workflow at $\varepsilon = 8$ for every geometry, "
    r"and, for the geometries whose solid mesh passed the standard OpenFOAM mesh checks, the mesh "
    r"size (median [interquartile range]), mesh quality (median [worst]) and the changes made to "
    r"the surface and the extrusion (median [range]). A point counts as moved when it is displaced "
    r"by more than \SI{0.02}{\milli\metre}; smoothing is the windowed-sinc smoothing of the surface, "
    r"and volume changes are of the lumen.}",
    r"\label{tab:meshing}",
    r"\begin{tabular}{@{}lccc@{}}",
    r"\toprule",
    r" & Human & Porcine & All \\",
    r"\midrule",
    *[line for line in table if line != ""],
    r"\bottomrule",
    r"\end{tabular}",
    r"\end{table}",
    "",
]))
(tables / "meshingNumbers.tex").write_text(
    "% Generated by scripts/meshing_results.py from data/meshing_ext8.csv; do not edit.\n"
    + "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in macros.items()))
print(f"Wrote tables/meshing.tex and {len(macros)} macros to tables/meshingNumbers.tex")
