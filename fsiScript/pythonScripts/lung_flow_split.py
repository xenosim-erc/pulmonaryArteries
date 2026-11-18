#!/usr/bin/env python3
"""Report the left/right lung flow split of finished run cases.

For each run case made by meshing.sh with lung groups, reads the flow rate
through the LPA and RPA outlets (postProcessing/fluid/{LPAFlow,RPAFlow}, all
start times merged, so restarted runs are handled) and, over the last complete
cardiac cycle, reports

- the mean flow through each lung and the left fraction of the total,
- the left fraction at peak systole (the time of the largest total outflow),
- the left fraction the outlet boundary conditions alone would give: the
  Windkessel totals are split by outlet area^(k/2), with k the flow-split
  exponent in windkessel/settings (2 if absent).

The cycle length is the "Heartbeats:" period in meshing.params.

Usage: lung_flow_split.py RUN_CASE [RUN_CASE ...] [--csv OUT.csv]
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np

trapezoid = getattr(np, "trapezoid", None) or np.trapz  # trapz was removed in NumPy 2


def read_flow(case: Path, name: str) -> np.ndarray:
    rows = {}
    files = sorted((case / "postProcessing" / "fluid" / name).glob("*/surfaceFieldValue*.dat"),
                   key=lambda f: float(f.parent.name))
    if not files:
        raise FileNotFoundError(f"no {name} output in {case}/postProcessing/fluid")
    for f in files:  # a later start time overrides the overlap of an earlier one
        for line in f.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            t, q = line.split()[:2]
            rows[float(t)] = float(q)
    data = np.array(sorted(rows.items()))
    return data


def cycle_period(case: Path) -> float:
    text = (case / "meshing.params").read_text()
    match = re.search(r"^Heartbeats:\s+\S+\s+x\s+([0-9.eE+-]+)\s+s", text, re.M)
    if not match:
        raise ValueError(f"no 'Heartbeats: N x T s' line in {case}/meshing.params")
    return float(match.group(1))


def boundary_condition_fraction(case: Path) -> float:
    exponent = 2.0
    settings = case / "windkessel" / "settings"
    if settings.exists():
        match = re.search(r"^flowSplitExponent=([0-9.eE+-]+)", settings.read_text(), re.M)
        if match:
            exponent = float(match.group(1))
    weight = {"L": 0.0, "R": 0.0}
    with (case / "lungSides.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            weight[row["side"]] += float(row["area"]) ** (exponent / 2)
    return weight["L"] / (weight["L"] + weight["R"])


def analyse(case: Path) -> dict:
    left, right = read_flow(case, "LPAFlow"), read_flow(case, "RPAFlow")
    times = np.intersect1d(left[:, 0], right[:, 0])
    ql = np.interp(times, left[:, 0], left[:, 1])
    qr = np.interp(times, right[:, 0], right[:, 1])
    period = cycle_period(case)
    cycles = int(np.floor(times[-1] / period + 1e-9))
    if cycles < 1:
        raise ValueError(f"{case}: only {times[-1]:.3g} s run, less than one cycle of {period} s")
    start, end = (cycles - 1) * period, cycles * period
    window = (times >= start - 1e-12) & (times <= end + 1e-12)
    t, ql, qr = times[window], ql[window], qr[window]
    mean_l = trapezoid(ql, t) / (t[-1] - t[0])
    mean_r = trapezoid(qr, t) / (t[-1] - t[0])
    peak = np.argmax(ql + qr)
    return dict(
        case=case.name,
        cycle=f"{start:.4g}-{end:.4g} s",
        mean_LPA_L_min=mean_l * 6e4,
        mean_RPA_L_min=mean_r * 6e4,
        left_fraction_mean=mean_l / (mean_l + mean_r),
        left_fraction_peak_systole=ql[peak] / (ql[peak] + qr[peak]),
        peak_systole_time_s=t[peak],
        left_fraction_boundary_conditions=boundary_condition_fraction(case),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("cases", nargs="+", type=Path)
    parser.add_argument("--csv", type=Path, help="write the results to this CSV")
    args = parser.parse_args()
    results, failed = [], 0
    for case in args.cases:
        try:
            results.append(analyse(case))
        except (OSError, ValueError, KeyError) as error:
            print(f"Skipping {case}: {error}", file=sys.stderr)
            failed += 1
    for r in results:
        print(f"{r['case']} ({r['cycle']}): LPA {r['mean_LPA_L_min']:.2f} L/min, RPA {r['mean_RPA_L_min']:.2f} L/min; "
              f"left {100 * r['left_fraction_mean']:.1f}% mean, {100 * r['left_fraction_peak_systole']:.1f}% at peak systole "
              f"(t = {r['peak_systole_time_s']:.3g} s); boundary conditions alone {100 * r['left_fraction_boundary_conditions']:.1f}%")
    if args.csv and results:
        with args.csv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(results[0]))
            writer.writeheader()
            for r in results:
                writer.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
        print(f"Wrote {args.csv}")
    return 1 if failed and not results else 0


if __name__ == "__main__":
    sys.exit(main())
