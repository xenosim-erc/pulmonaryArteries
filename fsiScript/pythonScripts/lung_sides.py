#!/usr/bin/env python3
"""Assign each outlet to the left or right lung.

Every centerline from pulmonary_centerlines.py runs from the inlet to one
outlet. For each pair of paths, the separation arc length is where one path
first moves more than half the local inscribed radius away from the other. The
outlets are clustered by single linkage on that length: the pair that stays
together longest is joined first, until two groups remain. The last two
groups part at the bifurcation of the main pulmonary artery. A branch that
leaves the main artery just before the bifurcation joins the side it stays
with longest.

A centerline with a straight jump longer than three local radii has left the
vessel (a broken Voronoi path; seen on 3 of the 75 cohort geometries), and
the topology cannot be trusted. The outlets are then split at the largest gap
between their centroids along --left-axis instead. This gap rule agreed with
the topology on all 72 cohort geometries with intact centerlines. When both
are available and disagree, the topology is used, with a warning.

The left lung is the group whose area-weighted outlet centroid lies further
along --left-axis. The default, +x, holds when the surfaces keep DICOM (LPS)
patient coordinates, in which +x is the patient's left; it was checked against
the anatomy of all the human and porcine cohort geometries that split cleanly.

A manual assignment (--override, a CSV of "name,side" with side L or R)
replaces the automatic one, for geometries that do not split cleanly. Every
outlet must be listed.

Writes a CSV of name, side (L or R) and area, and prints the split of outlet
area between the lungs: with the Windkessel totals divided by area, this is the
flow split the boundary conditions alone would give (exponent 2).

Exit status: 0 assigned (check any warning), 3 not assigned (fewer than two
outlets; the case is still usable without lung groups), 2 bad input.
"""

from __future__ import annotations

import argparse
import csv
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

AXES = {"+x": (0, 1), "-x": (0, -1), "+y": (1, 1), "-y": (1, -1), "+z": (2, 1), "-z": (2, -1)}
JUMP_RADII = 3.0  # a centerline step longer than this many radii is a break


def read_centerlines(path: Path) -> dict[int, np.ndarray]:
    paths: dict[int, list] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            paths.setdefault(int(row["branch"]), []).append(
                [float(row[k]) for k in ("x", "y", "z", "distance", "radius")])
    return {b: np.asarray(p) for b, p in paths.items()}


def read_profiles(path: Path):
    outlets, inlet = {}, None
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            centre = np.array([float(row[k]) for k in "xyz"])
            if row["name"] == "inlet":
                inlet = centre
            else:
                outlets[row["name"]] = (float(row["area"]), centre)
    if inlet is None or not outlets:
        raise ValueError(f"{path} has no inlet or no outlets")
    return outlets, inlet


def separation(a: np.ndarray, b: np.ndarray) -> float:
    d = np.sqrt(((a[:, None, :3] - b[None, :, :3]) ** 2).sum(-1)).min(axis=1)
    away = np.nonzero(d > 0.5 * a[:, 4])[0]
    return a[away[0], 3] if away.size else a[-1, 3]


def broken(paths) -> bool:
    for p in paths.values():
        step = np.linalg.norm(np.diff(p[:, :3], axis=0), axis=1)
        if (step > JUMP_RADII * np.maximum(p[:-1, 4], p[1:, 4])).any():
            return True
    return False


def topology_groups(paths, path_outlet):
    """Two outlet groups by single linkage on the separation arc length."""
    branches = sorted(paths)
    parent = {b: b for b in branches}

    def find(b):
        while parent[b] != b:
            parent[b] = parent[parent[b]]
            b = parent[b]
        return b

    pairs = sorted(((min(separation(paths[i], paths[j]), separation(paths[j], paths[i])), i, j)
                    for i, j in combinations(branches, 2)), reverse=True)
    clusters, split = len(branches), None
    for together, i, j in pairs:
        a, b = find(i), find(j)
        if a == b:
            continue
        if clusters == 2:
            split = together
            break
        parent[a] = b
        clusters -= 1
    groups = {}
    for b in branches:
        groups.setdefault(find(b), []).append(path_outlet[b])
    return list(groups.values()), split


def gap_groups(outlets, axis):
    """Two outlet groups either side of the largest gap along the axis."""
    names = sorted(outlets, key=lambda n: outlets[n][1][axis])
    positions = np.array([outlets[n][1][axis] for n in names])
    k = int(np.diff(positions).argmax())
    return [names[:k + 1], names[k + 1:]]


def automatic_sides(paths, outlets, left_axis):
    """Return ({name: 'L'|'R'}, message) or (None, reason)."""
    axis, sign = AXES[left_axis]
    if len(outlets) < 2:
        return None, "fewer than two outlets"

    def centroid(group):
        areas = np.array([outlets[n][0] for n in group])
        return (areas[:, None] * np.array([outlets[n][1] for n in group])).sum(0) / areas.sum()

    def label(groups):
        left = 0 if sign * centroid(groups[0])[axis] > sign * centroid(groups[1])[axis] else 1
        return {n: ("L" if k == left else "R") for k, g in enumerate(groups) for n in g}

    by_gap = label(gap_groups(outlets, axis))
    path_outlet = {b: min(outlets, key=lambda n: np.linalg.norm(outlets[n][1] - p[-1, :3]))
                   for b, p in paths.items()}
    if sorted(path_outlet.values()) != sorted(outlets):
        return by_gap, ("the centerline ends do not match the outlets one-to-one; "
                        "split at the largest gap between outlets instead (check it)")
    if broken(paths):
        return by_gap, ("a centerline jumps by more than "
                        f"{JUMP_RADII:g} radii (broken); split at the largest gap "
                        "between outlets instead (check it)")
    groups, split = topology_groups(paths, path_outlet)
    by_topology = label(groups)
    message = f"main pulmonary artery bifurcation {split:.1f} mm along the centerlines"
    if by_topology != by_gap:
        message += "; WARNING: the largest-gap split disagrees, check this case"
    return by_topology, message


def read_override(path: Path, outlets):
    sides = {}
    with path.open(newline="") as handle:
        for row in csv.reader(handle):
            if not row or row[0].strip().startswith("#") or row[0].strip() == "name":
                continue
            if len(row) < 2:
                raise ValueError(f"{path}: expected 'name,side', got {','.join(row)!r}")
            name, side = row[0].strip(), row[1].strip().upper()
            if side not in ("L", "R"):
                raise ValueError(f"{path}: side of {name} must be L or R, not {row[1]!r}")
            sides[name] = side
    if sorted(sides) != sorted(outlets):
        missing = sorted(set(outlets) - set(sides))
        extra = sorted(set(sides) - set(outlets))
        raise ValueError(f"{path} must list every outlet once (missing {missing}, unknown {extra})")
    return sides


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("centerlines_csv", type=Path)
    parser.add_argument("profiles_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--left-axis", choices=sorted(AXES), default="+x",
                        help="direction of the patient's left in the surface coordinates (default +x, LPS)")
    parser.add_argument("--override", type=Path, help="CSV of name,side (L or R) for every outlet")
    args = parser.parse_args()
    try:
        outlets, _ = read_profiles(args.profiles_csv)
        if args.override:
            sides, how = read_override(args.override, outlets), f"from {args.override}"
        else:
            sides, how = automatic_sides(read_centerlines(args.centerlines_csv), outlets, args.left_axis)
            if sides is None:
                print(f"Warning: lung sides not assigned: {how}", file=sys.stderr)
                return 3
            if "check" in how:
                print(f"Warning: lung sides: {how}", file=sys.stderr)
            how = f"automatic, left = {args.left_axis}; {how}"
    except (OSError, ValueError, KeyError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2

    with args.output_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["name", "side", "area"])
        for name in sorted(outlets, key=lambda n: (sides[n], int("".join(filter(str.isdigit, n)) or 0))):
            writer.writerow([name, sides[name], f"{outlets[name][0]:.6g}"])
    area = {s: sum(outlets[n][0] for n in outlets if sides[n] == s) for s in "LR"}
    total = area["L"] + area["R"]
    for s, label in (("L", "Left"), ("R", "Right")):
        names = sorted((n for n in outlets if sides[n] == s), key=lambda n: int("".join(filter(str.isdigit, n)) or 0))
        print(f"{label} lung: {len(names)} outlet(s), {100 * area[s] / total:.1f}% of outlet area: {' '.join(names)}")
    print(f"Lung sides {how}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
