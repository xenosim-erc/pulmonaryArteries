#!/usr/bin/env python3
"""Refine the cfMesh dictionary around outlets too small for the cell size.

An outlet only a few cells across is meshed shut: the opening never appears in
the volume mesh, so the case silently loses a boundary condition. Each profile
whose radius is small relative to maxCellSize therefore gets a refinement
sphere, sized from its own radius rather than placed by hand, so the mesh
dictionary stays independent of any particular anatomy.

Refinement is deliberately shallow. A patch of wall meshed much finer than its
surroundings leaves a size transition there, and cfMesh's staircase across that
transition is sharp enough to damage the solid extrusion that follows, so
--maximum-refinement-levels caps how far below maxCellSize this may go.
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


def build_block(refined, radius_factor: float) -> str:
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
            "refinement sphere radius as a multiple of the profile radius "
            "(default: 4). The sphere is centred on the extended rim, and a "
            "flow extension of one diameter is itself two radii long, so the "
            "sphere must span at least that to reach the cap's feeder vessel. "
            "Shorten the extensions and this should come down with them"
        ),
    )
    parser.add_argument(
        "--maximum-refinement-levels",
        type=int,
        default=1,
        help=(
            "most octree levels any refinement may add below maxCellSize "
            "(default: 1, so refined cells are at most half the global size)"
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


    # A profile only needs refining when the cell size it wants is finer than
    # the global one; larger vessels are already resolved.
    #
    # The requested size is snapped to one of cfMesh's octree levels rather than
    # passed through. cfMesh can only halve, so it satisfies a request by taking
    # the first level at or below it: asking for r/3 = 0.4997 against a global
    # size of 1 does not give 0.5, it gives 0.25, which is twice the intended
    # resolution and eight times the cells. Levels are then capped, because a
    # refined patch of wall that is much finer than its surroundings leaves a
    # size transition sharp enough to behave as a separate surface.
    refined: list[tuple[dict[str, float], float]] = []
    for profile in profiles:
        wanted = profile["radius"] / args.cells_per_radius
        if wanted >= max_cell_size:
            continue
        levels = math.ceil(math.log2(max_cell_size / wanted))
        levels = min(max(levels, 1), args.maximum_refinement_levels)
        cell_size = max(max_cell_size / (2.0**levels), args.minimum_cell_size)
        refined.append((profile, cell_size))

    if not refined:
        updated = remove_existing_block(text)
        if updated != text:
            write_atomically(mesh_dict, updated)
        print(
            f"Nothing is narrow relative to maxCellSize {max_cell_size:g}; "
            "no refinement added"
        )
        return 0

    updated = remove_existing_block(text).rstrip("\n") + "\n\n"
    updated += build_block(refined, args.sphere_radius_factor)
    write_atomically(mesh_dict, updated)

    sizes = [cell_size for _, cell_size in refined]
    print(
        f"Refined {len(refined)} small outlet(s) of {len(profiles)} to cell "
        f"sizes {min(sizes):.3g}-{max(sizes):.3g}, against maxCellSize "
        f"{max_cell_size:g}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, OSError, RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
