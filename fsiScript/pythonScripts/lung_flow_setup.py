#!/usr/bin/env python3
"""Group the outlet patches by lung and write flow and pressure function objects.

Reads the lung sides from lung_sides.py (name, side, area) and

1. adds every outlet patch of the fluid boundary file to the patch group LPA
   or RPA (an "inGroups" entry), and
2. writes a dictionary of function objects, included by the case's
   system/controlDict, that record at every time step the volumetric flow rate
   (sum of phi, <name>Flow) and the area-averaged kinematic pressure (p, in
   m^2/s^2; <name>Pressure) through the inlet, each lung's outlets and each
   outlet.

surfaceFieldValue matches patch groups by name, so "names (LPA)" sums over all
the left-lung outlets. Output: postProcessing/fluid/<name>/<startTime>/
surfaceFieldValue.dat. The outlet rings of the wall are fixed, so phi there is
the flow rate of the blood. The inlet flow is negative (into the domain), so
inletFlow + LPAFlow + RPAFlow is the mass-balance error.

Usage: lung_flow_setup.py SIDES_CSV FLUID_BOUNDARY_FILE FUNCTIONS_FILE
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

GROUP = {"L": "LPA", "R": "RPA"}


def read_sides(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        return {row["name"]: row["side"] for row in csv.DictReader(handle)}


def add_groups(boundary: Path, sides: dict[str, str]) -> None:
    text = boundary.read_text()
    for name, side in sides.items():
        block = re.search(rf"(\n[ \t]*{re.escape(name)}[ \t]*\n[ \t]*\{{)(.*?)(\n[ \t]*\}})", text, re.S)
        if block is None:
            raise ValueError(f"patch {name} not found in {boundary}")
        body = block.group(2)
        group = GROUP[side]
        existing = re.search(r"\n([ \t]*)inGroups[ \t]+(?:List<word>[ \t]*)?\d*\((.*?)\);", body)
        if existing:
            groups = existing.group(2).split()
            if group not in groups:
                groups.append(group)
            new = f"\n{existing.group(1)}inGroups        List<word> {len(groups)}({' '.join(groups)});"
            body = body[:existing.start()] + new + body[existing.end():]
        else:
            indent = re.search(r"\n([ \t]*)type", body).group(1)
            body = re.sub(r"(\n[ \t]*type[^\n]*;)", rf"\1\n{indent}inGroups        List<word> 1({group});", body, count=1)
        text = text[:block.start(2)] + body + text[block.end(2):]
    boundary.write_text(text)


def function_object(name: str, selection: str, operation: str = "sum", field: str = "phi") -> str:
    return f"""{name}
{{
    type            surfaceFieldValue;
    libs            (fieldFunctionObjects);
    region          fluid;
    regionType      patch;
    names           ({selection});
    operation       {operation};
    fields          ({field});
    writeFields     false;
    writeControl    timeStep;
    writeInterval   1;
    log             false;
}}
"""


def write_functions(path: Path, sides: dict[str, str]) -> None:
    order = sorted(sides, key=lambda n: int("".join(filter(str.isdigit, n)) or 0))
    parts = [
        "// Flow rate and mean pressure through the inlet, each lung's outlets and each\n"
        "// outlet. Written by meshing.sh (pythonScripts/lung_flow_setup.py) and\n"
        "// included by system/controlDict.\n"
        "// " + "; ".join(f"{GROUP[s]}: {' '.join(n for n in order if sides[n] == s)}" for s in "LR") + "\n",
    ]
    for name in ["inlet", "LPA", "RPA", *order]:
        parts.append(function_object(f"{name}Flow", name))
        parts.append(function_object(f"{name}Pressure", name, "areaAverage", "p"))
    path.write_text("\n".join(parts))


def main() -> int:
    if len(sys.argv) != 4:
        print(__doc__, file=sys.stderr)
        return 2
    sides_csv, boundary, functions = map(Path, sys.argv[1:])
    try:
        sides = read_sides(sides_csv)
        add_groups(boundary, sides)
        write_functions(functions, sides)
    except (OSError, ValueError, KeyError, AttributeError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    n = {s: sum(1 for v in sides.values() if v == s) for s in "LR"}
    print(f"Grouped {n['L']} outlet(s) as LPA and {n['R']} as RPA in {boundary}")
    print(f"Wrote inlet, lung and outlet flow and pressure function objects to {functions}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
