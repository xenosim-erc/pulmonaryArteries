#!/usr/bin/env python3
"""Plot the inlet pressure-volume curves and power measures across cases.

Reads fsiScript/results/summary.csv and each results/<case>/pv_inlet.csv
(written by extract_results.py) and writes, in results/plots/:

  pv_ejection.png   inlet pressure against volume ejected so far in the beat:
                    the ejection limb of a right-ventricular PV loop; the area
                    under each curve is the hydraulic work per beat
  pv_lumen.png      inlet pressure against the change in lumen volume of the
                    3D arteries: a closed loop whose slope is their compliance
  power.png         oscillatory power fraction of each case, FSI against its
                    outlet-Windkessel-only baseline

Cases are coloured, and given a line style, by species (the "animal" column).
--normalise divides volumes by each case's stroke volume, so cases of
different size can be overlaid. Runs with the system python3 (needs
matplotlib and numpy only).

Usage: plot_results.py [--results-dir DIR] [--cases NAME ...] [--normalise]
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

FSI_DIR = Path(__file__).resolve().parent.parent
# Categorical slots in fixed order (validated for colour-vision deficiency);
# each species also gets its own line style, so colour is never the only cue.
SPECIES = {"human": ("#2a78d6", "-", "o"), "pig": ("#eb6834", "--", "s")}
OTHER = ("#1baf7a", ":", "^")
INK, MUTED, GRID = "#1f1f1e", "#6b6a63", "#e4e3dd"


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def number(value) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def style(ax, title, xlabel, ylabel):
    ax.set_title(title, loc="left", color=INK, fontsize=11)
    ax.set_xlabel(xlabel, color=MUTED)
    ax.set_ylabel(ylabel, color=MUTED)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED)


def species_legend(ax, present):
    handles = [plt.Line2D([], [], color=SPECIES.get(s, OTHER)[0], linestyle=SPECIES.get(s, OTHER)[1],
                          linewidth=2, label=s) for s in present]
    ax.legend(handles=handles, frameon=False, labelcolor=INK)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=__doc__.split("\n\n", 1)[1])
    parser.add_argument("--results-dir", type=Path, default=FSI_DIR / "results")
    parser.add_argument("--cases", nargs="+", help="only these cases (default: every row of summary.csv)")
    parser.add_argument("--normalise", action="store_true", help="volumes as a fraction of stroke volume")
    args = parser.parse_args()

    summary_path = args.results_dir / "summary.csv"
    if not summary_path.is_file():
        print(f"error: no {summary_path}; run extract_results.py first", file=sys.stderr)
        return 1
    rows = [r for r in read_csv(summary_path) if not args.cases or r["case"] in args.cases]
    plots = args.results_dir / "plots"
    plots.mkdir(exist_ok=True)

    fig_e, ax_e = plt.subplots(figsize=(7, 4.5), dpi=150)
    fig_l, ax_l = plt.subplots(figsize=(7, 4.5), dpi=150)
    present, plotted = [], 0
    for row in rows:
        pv_path = args.results_dir / row["case"] / "pv_inlet.csv"
        if not pv_path.is_file():
            continue
        pv = read_csv(pv_path)
        p = [number(r["p_inlet_mmHg"]) for r in pv]
        ejected = [number(r["V_ejected_mL"]) for r in pv]
        lumen = [number(r["dVolume_mL"]) for r in pv]
        scale = 1.0
        if args.normalise:
            stroke = ejected[-1] if ejected and ejected[-1] > 0 else math.nan
            scale = 1.0 / stroke
        species = (row.get("animal") or "other").lower()
        color, line, _ = SPECIES.get(species, OTHER)
        if species not in present:
            present.append(species)
        ax_e.plot([v * scale for v in ejected], p, color=color, linestyle=line, linewidth=2)
        if not all(math.isnan(v) for v in lumen):
            ax_l.plot([v * scale for v in lumen], p, color=color, linestyle=line, linewidth=2)
        # Direct label at the end of each ejection curve.
        ax_e.annotate(row["case"], (ejected[-1] * scale, p[-1]), xytext=(4, 0), textcoords="offset points",
                      fontsize=7, color=MUTED, va="center")
        plotted += 1
    if not plotted:
        print("error: no pv_inlet.csv found for the selected cases", file=sys.stderr)
        return 1

    unit = "fraction of stroke volume" if args.normalise else "mL"
    style(ax_e, "Inlet pressure against ejected volume (last cycle)", f"Volume ejected ({unit})",
          "Inlet pressure (mmHg)")
    style(ax_l, "Inlet pressure against lumen volume of the 3D arteries", f"Change in lumen volume ({unit})",
          "Inlet pressure (mmHg)")
    for ax in (ax_e, ax_l):
        species_legend(ax, present)
    fig_e.tight_layout()
    fig_l.tight_layout()
    fig_e.savefig(plots / "pv_ejection.png")
    fig_l.savefig(plots / "pv_lumen.png")

    # Oscillatory power fraction: FSI (filled) against its Windkessel-only
    # baseline (hollow), one row per case, grouped by species.
    ranked = sorted((r for r in rows if not math.isnan(number(r.get("oscillatoryPowerFraction")))),
                    key=lambda r: ((r.get("animal") or "other").lower(), r["case"]))
    if ranked:
        fig, ax = plt.subplots(figsize=(7, 0.45 * len(ranked) + 1.4), dpi=150)
        for i, r in enumerate(ranked):
            color, _, marker = SPECIES.get((r.get("animal") or "other").lower(), OTHER)
            fsi, base = number(r["oscillatoryPowerFraction"]), number(r.get("baseline_oscillatoryPowerFraction"))
            if not math.isnan(base):
                ax.plot([base, fsi], [i, i], color=GRID, linewidth=2, zorder=1)
                ax.scatter([base], [i], s=64, facecolor="white", edgecolor=color, linewidth=2, marker=marker,
                           zorder=2)
            ax.scatter([fsi], [i], s=64, color=color, edgecolor="white", linewidth=2, marker=marker, zorder=3)
        ax.set_yticks(range(len(ranked)), [r["case"] for r in ranked])
        style(ax, "Oscillatory share of inlet hydraulic power", "Oscillatory power / total power", "")
        ax.set_xlim(left=0)
        handles = [plt.Line2D([], [], color=INK, marker="o", linestyle="", markersize=7, label="FSI"),
                   plt.Line2D([], [], color=INK, marker="o", markerfacecolor="white", linestyle="", markersize=7,
                              label="outlet Windkessels only")]
        handles += [plt.Line2D([], [], color=SPECIES.get(s, OTHER)[0], marker=SPECIES.get(s, OTHER)[2],
                               linestyle="", markersize=7, label=s) for s in present]
        ax.legend(handles=handles, frameon=False, labelcolor=INK, loc="best", fontsize=8)
        fig.tight_layout()
        fig.savefig(plots / "power.png")
    print(f"wrote {plots} ({plotted} case(s))")
    return 0


if __name__ == "__main__":
    sys.exit(main())
