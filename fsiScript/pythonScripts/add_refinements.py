#!/usr/bin/env python3
"""Refine the cfMesh dictionary wherever the vessel is narrow for the cell size.

Two things need resolving beyond the global maxCellSize.

An outlet cap must be resolved by enough cells for autoPatch to separate it from
the vessel wall. A profile whose radius is comparable to maxCellSize is only a
few cells across, its cap merges into the wall, and the merged faces then have
no counterpart on the uncapped surface that supplies the wall thickness. Each
such profile gets a sphere.

A narrow branch anywhere along the tree has the same problem in a different
form: a Cartesian mesh snapped to a tube only six or seven cells across leaves
faceted corners on the wall, and offsetting those corners along their point
normals produces warped and incorrectly oriented faces in the extruded solid.
Centerlines carry a radius at every point, so the branches are refined by cone
segments following the vessel itself.

Both are derived from the geometry rather than placed by hand, so the mesh
dictionary stays independent of any particular anatomy.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import stat
import sys
import tempfile
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


MAX_CELL_SIZE = re.compile(r"(?m)^[ \t]*maxCellSize[ \t]+(?P<size>[0-9.eE+-]+)[ \t]*;")


def positive_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise argparse.ArgumentTypeError("must be finite and greater than 0")
    return number


def read_max_cell_size(text: str) -> float:
    match = MAX_CELL_SIZE.search(text)
    if match is None:
        raise RuntimeError("could not find maxCellSize in the mesh dictionary")
    size = float(match.group("size"))
    if not math.isfinite(size) or size <= 0.0:
        raise RuntimeError(f"maxCellSize is not a positive number: {size}")
    return size


def read_profiles(filename: Path) -> list[dict[str, float]]:
    with filename.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"index", "radius", "x", "y", "z"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise RuntimeError(
                f"profile CSV is missing columns: {', '.join(sorted(missing))}"
            )
        profiles = []
        for row_number, row in enumerate(reader, start=2):
            try:
                profile = {
                    "index": int(row["index"]),
                    "radius": float(row["radius"]),
                    "centre": [float(row[axis]) for axis in ("x", "y", "z")],
                }
            except (TypeError, ValueError) as exc:
                raise RuntimeError(
                    f"invalid numeric data on profile CSV row {row_number}"
                ) from exc
            if profile["radius"] <= 0.0:
                raise RuntimeError(
                    f"non-positive profile radius on row {row_number}"
                )
            profiles.append(profile)
    if not profiles:
        raise RuntimeError(f"no profiles were found in {filename}")
    return profiles


def read_centerlines(filename: Path):
    """Centerline points, their vessel radius, and the branch each belongs to."""
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(filename))
    reader.Update()
    centerlines = reader.GetOutput()
    radius = centerlines.GetPointData().GetArray("MaximumInscribedSphereRadius")
    if radius is None:
        raise RuntimeError(
            "centerlines are missing the MaximumInscribedSphereRadius array"
        )

    branches = []
    lines = centerlines.GetLines()
    lines.InitTraversal()
    identifiers = vtk.vtkIdList()
    while lines.GetNextCell(identifiers):
        ids = [identifiers.GetId(i) for i in range(identifiers.GetNumberOfIds())]
        points = np.array([centerlines.GetPoint(i) for i in ids])
        radii = np.array([radius.GetTuple1(i) for i in ids])
        if len(points) >= 2:
            branches.append((points, radii))
    if not branches:
        raise RuntimeError(f"no centerline branches were found in {filename}")
    return branches


def cone_refinements(branches, max_cell_size, cells_per_radius, radius_factor):
    """One cone per run of centerline points that wants the same cell size.

    Required cell sizes are quantised to the octree levels cfMesh works in, so a
    branch becomes a handful of cones instead of one sphere per centerline point.
    VMTK traces every branch from the same source, so the shared trunk is
    repeated in each line; identical segments are emitted only once.
    """
    cones = []
    seen = set()
    for points, radii in branches:
        wanted = radii / cells_per_radius
        # Octree level: 0 means the global size is already fine enough.
        levels = np.where(
            wanted < max_cell_size,
            np.ceil(np.log2(max_cell_size / np.maximum(wanted, 1e-12))),
            0,
        ).astype(int)

        start = 0
        for index in range(1, len(levels) + 1):
            if index < len(levels) and levels[index] == levels[start]:
                continue
            level = int(levels[start])
            if level > 0:
                segment = slice(start, index)
                first, last = points[start], points[index - 1]
                if not np.allclose(first, last):
                    key = (
                        level,
                        tuple(np.round(first, 3)),
                        tuple(np.round(last, 3)),
                    )
                    if key not in seen:
                        seen.add(key)
                        cones.append(
                            {
                                "p0": first,
                                "p1": last,
                                "radius0": radius_factor * radii[start],
                                "radius1": radius_factor * radii[index - 1],
                                "cellSize": max_cell_size / (2.0**level),
                            }
                        )
            start = index
    return cones


def remove_existing_block(text: str) -> str:
    """Strip any objectRefinements block so this script owns the whole entry."""
    match = re.search(r"(?m)^[ \t]*objectRefinements\b", text)
    if match is None:
        return text
    opening = text.find("{", match.end())
    if opening < 0:
        raise RuntimeError("objectRefinements is present but has no opening brace")
    depth = 0
    for position in range(opening, len(text)):
        if text[position] == "{":
            depth += 1
        elif text[position] == "}":
            depth -= 1
            if depth == 0:
                return text[: match.start()] + text[position + 1 :].lstrip("\n")
    raise RuntimeError("objectRefinements block is not brace balanced")


def build_block(refined, radius_factor: float, cones) -> str:
    lines = ["objectRefinements", "{"]
    for profile, cell_size in refined:
        centre = " ".join(f"{value:.6g}" for value in profile["centre"])
        lines += [
            f"    smallOutlet{profile['index']}",
            "    {",
            "        type        sphere;",
            f"        cellSize    {cell_size:.6g};",
            f"        centre      ({centre});",
            f"        radius      {radius_factor * profile['radius']:.6g};",
            "    }",
        ]
    for index, cone in enumerate(cones):
        p0 = " ".join(f"{value:.6g}" for value in cone["p0"])
        p1 = " ".join(f"{value:.6g}" for value in cone["p1"])
        lines += [
            f"    narrowVessel{index}",
            "    {",
            "        type        cone;",
            f"        cellSize    {cone['cellSize']:.6g};",
            f"        p0          ({p0});",
            f"        radius0     {cone['radius0']:.6g};",
            f"        p1          ({p1});",
            f"        radius1     {cone['radius1']:.6g};",
            "    }",
        ]
    lines += ["}", ""]
    return "\n".join(lines)


def write_atomically(path: Path, text: str) -> None:
    original_mode = stat.S_IMODE(path.stat().st_mode)
    temporary_name = ""
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            delete=False,
        ) as temporary:
            temporary.write(text)
            temporary.flush()
            os.fsync(temporary.fileno())
            temporary_name = temporary.name
        os.chmod(temporary_name, original_mode)
        os.replace(temporary_name, path)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh_dict", type=Path)
    parser.add_argument("profiles_csv", type=Path)
    parser.add_argument(
        "--centerlines",
        type=Path,
        default=None,
        help=(
            "centerline VTP enabling narrow-vessel cone refinement. OFF by "
            "default: on the supplied geometries the refinement-level "
            "transitions it introduces on the wall patch damage the extruded "
            "solid far more than the faceting it removes"
        ),
    )
    parser.add_argument(
        "--cells-per-radius",
        type=positive_float,
        default=3.0,
        help="cells wanted across a vessel radius (default: 3)",
    )
    parser.add_argument(
        "--sphere-radius-factor",
        type=positive_float,
        default=4.0,
        help=(
            "refinement sphere radius as a multiple of the profile radius. The "
            "sphere is centred on the extended rim, and a flow extension is "
            "itself one diameter long, so a factor of 2 would reach only back "
            "to the original rim and cover none of the feeder vessel"
        ),
    )
    parser.add_argument(
        "--cone-radius-factor",
        type=positive_float,
        default=1.5,
        help=(
            "narrow-vessel cone radius as a multiple of the local vessel radius "
            "(default: 1.5). Only the wall needs covering, so this stays much "
            "tighter than the outlet spheres, which must also span the extension"
        ),
    )
    parser.add_argument(
        "--minimum-cell-size",
        type=positive_float,
        default=0.0,
        help="never refine below this cell size; 0 means no floor",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    mesh_dict = args.mesh_dict.expanduser().resolve()
    profiles_csv = args.profiles_csv.expanduser().resolve()
    for filename in (mesh_dict, profiles_csv):
        if not filename.is_file():
            raise FileNotFoundError(filename)

    text = mesh_dict.read_text(encoding="utf-8")
    max_cell_size = read_max_cell_size(text)
    profiles = read_profiles(profiles_csv)

    cones = []
    if args.centerlines is not None:
        centerlines = args.centerlines.expanduser().resolve()
        if not centerlines.is_file():
            raise FileNotFoundError(centerlines)
        cones = cone_refinements(
            read_centerlines(centerlines),
            max_cell_size,
            args.cells_per_radius,
            args.cone_radius_factor,
        )

    # A profile only needs refining when the cell size it wants is finer than
    # the global one; larger vessels are already resolved.
    refined: list[tuple[dict[str, float], float]] = []
    for profile in profiles:
        wanted = profile["radius"] / args.cells_per_radius
        if wanted >= max_cell_size:
            continue
        refined.append((profile, max(wanted, args.minimum_cell_size)))

    if not refined and not cones:
        updated = remove_existing_block(text)
        if updated != text:
            write_atomically(mesh_dict, updated)
        print(
            f"Nothing is narrow relative to maxCellSize {max_cell_size:g}; "
            "no refinement added"
        )
        return 0

    updated = remove_existing_block(text).rstrip("\n") + "\n\n"
    updated += build_block(refined, args.sphere_radius_factor, cones)
    write_atomically(mesh_dict, updated)

    sizes = [cell_size for _, cell_size in refined] + [c["cellSize"] for c in cones]
    print(
        f"Refined {len(refined)} small outlet(s) of {len(profiles)} and "
        f"{len(cones)} narrow vessel segment(s) to cell sizes "
        f"{min(sizes):.3g}-{max(sizes):.3g}, against maxCellSize {max_cell_size:g}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, OSError, RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
