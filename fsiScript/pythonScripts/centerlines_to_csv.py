#!/usr/bin/env python3
"""Export VMTK centerline points and radii from VTP to CSV."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import vtk


RADIUS_ARRAY_NAME = "MaximumInscribedSphereRadius"


def read_centerlines(filename: Path) -> vtk.vtkPolyData:
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(filename))
    reader.Update()
    centerlines = reader.GetOutput()
    if centerlines.GetNumberOfLines() == 0:
        raise RuntimeError(f"No centerline branches were found in {filename}")
    return centerlines


def write_csv(centerlines: vtk.vtkPolyData, filename: Path) -> None:
    radius_array = centerlines.GetPointData().GetArray(RADIUS_ARRAY_NAME)
    if radius_array is None:
        raise RuntimeError(
            f"Required point-data array '{RADIUS_ARRAY_NAME}' was not found"
        )

    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(
            ["branch", "point", "x", "y", "z", "distance", "radius", "diameter"]
        )

        # Traverse each vtkPolyLine once. Point IDs refer back to coordinates
        # and MaximumInscribedSphereRadius values in the containing PolyData.
        lines = centerlines.GetLines()
        lines.InitTraversal()
        point_ids = vtk.vtkIdList()
        branch = 0
        while lines.GetNextCell(point_ids):
            distance = 0.0
            previous_point: tuple[float, float, float] | None = None

            for point_index in range(point_ids.GetNumberOfIds()):
                point_id = point_ids.GetId(point_index)
                point = centerlines.GetPoint(point_id)
                radius = radius_array.GetTuple1(point_id)

                if previous_point is not None:
                    distance += math.dist(previous_point, point)

                writer.writerow(
                    [
                        branch,
                        point_index,
                        point[0],
                        point[1],
                        point[2],
                        distance,
                        radius,
                        2.0 * radius,
                    ]
                )
                previous_point = point
            branch += 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_vtp", type=Path, help="VMTK centerlines VTP file")
    parser.add_argument("output_csv", type=Path, help="destination CSV file")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_path = args.input_vtp.expanduser().resolve()
    output_path = args.output_csv.expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    write_csv(read_centerlines(input_path), output_path)
    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, OSError, RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
