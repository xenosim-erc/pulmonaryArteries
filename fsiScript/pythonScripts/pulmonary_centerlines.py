#!/usr/bin/env python3
"""Convert an open pulmonary-artery STL to VTP and extract centerlines.

The largest open boundary (by planar profile area) is used as the source;
every other open boundary is used as a target.  VTK performs all surface I/O
and preparation.  VMTK's Python bindings perform the centerline extraction.
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict, deque
from pathlib import Path
from typing import Sequence

import vtk


Point = tuple[float, float, float]


def read_stl(filename: Path) -> vtk.vtkPolyData:
    reader = vtk.vtkSTLReader()
    reader.SetFileName(str(filename))
    reader.Update()
    if reader.GetOutput().GetNumberOfCells() == 0:
        raise RuntimeError(f"No triangles could be read from {filename}")
    return reader.GetOutput()


def prepare_surface(surface: vtk.vtkPolyData) -> vtk.vtkPolyData:
    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(surface)
    clean.PointMergingOn()
    clean.Update()

    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputConnection(clean.GetOutputPort())
    triangles.PassLinesOff()
    triangles.PassVertsOff()
    triangles.Update()

    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(triangles.GetOutputPort())
    normals.ConsistencyOn()
    normals.AutoOrientNormalsOn()
    normals.SplittingOff()
    normals.Update()

    output = vtk.vtkPolyData()
    output.DeepCopy(normals.GetOutput())
    return output


def boundary_loops(surface: vtk.vtkPolyData) -> list[list[Point]]:
    """Return ordered point loops formed by the surface boundary edges."""
    edges = vtk.vtkFeatureEdges()
    edges.SetInputData(surface)
    edges.BoundaryEdgesOn()
    edges.FeatureEdgesOff()
    edges.NonManifoldEdgesOff()
    edges.ManifoldEdgesOff()
    edges.Update()
    boundary = edges.GetOutput()

    adjacency: dict[int, set[int]] = defaultdict(set)
    for cell_id in range(boundary.GetNumberOfCells()):
        cell = boundary.GetCell(cell_id)
        if cell.GetNumberOfPoints() != 2:
            continue
        a, b = cell.GetPointId(0), cell.GetPointId(1)
        adjacency[a].add(b)
        adjacency[b].add(a)

    loops: list[list[Point]] = []
    unseen = set(adjacency)
    while unseen:
        root = next(iter(unseen))
        component: set[int] = set()
        queue = deque([root])
        while queue:
            point_id = queue.popleft()
            if point_id in component:
                continue
            component.add(point_id)
            queue.extend(adjacency[point_id] - component)
        unseen -= component

        bad = [point_id for point_id in component if len(adjacency[point_id]) != 2]
        if bad:
            raise RuntimeError(
                "A boundary is not a simple closed loop. The STL may contain "
                "non-manifold edges or cracks."
            )
        ordered = [root]
        previous = -1
        current = root
        while True:
            following = next(p for p in adjacency[current] if p != previous)
            if following == root:
                break
            ordered.append(following)
            previous, current = current, following
            if len(ordered) > len(component):
                raise RuntimeError("Failed to order an open boundary.")
        loops.append([tuple(boundary.GetPoint(i)) for i in ordered])
    return loops


def profile_area_and_centroid(points: Sequence[Point]) -> tuple[float, Point]:
    """Return projected polygon area and vertex centroid for a near-planar rim."""
    if len(points) < 3:
        return 0.0, (0.0, 0.0, 0.0)
    nx = ny = nz = 0.0
    for p, q in zip(points, (*points[1:], points[0])):
        nx += (p[1] - q[1]) * (p[2] + q[2])
        ny += (p[2] - q[2]) * (p[0] + q[0])
        nz += (p[0] - q[0]) * (p[1] + q[1])
    area = 0.5 * math.sqrt(nx * nx + ny * ny + nz * nz)
    count = float(len(points))
    centroid = tuple(sum(p[axis] for p in points) / count for axis in range(3))
    return area, centroid  # type: ignore[return-value]


def cap_surface_vtk(surface: vtk.vtkPolyData) -> vtk.vtkPolyData:
    fill = vtk.vtkFillHolesFilter()
    fill.SetInputData(surface)
    fill.SetHoleSize(sys.float_info.max)
    fill.Update()
    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputConnection(fill.GetOutputPort())
    triangles.Update()
    output = vtk.vtkPolyData()
    output.DeepCopy(triangles.GetOutput())
    return output


def write_vtp(surface: vtk.vtkPolyData, filename: Path) -> None:
    filename.parent.mkdir(parents=True, exist_ok=True)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(filename))
    writer.SetInputData(surface)
    writer.SetDataModeToBinary()
    if writer.Write() != 1:
        raise RuntimeError(f"Could not write {filename}")


def write_stl(surface: vtk.vtkPolyData, filename: Path) -> None:
    """Write a triangulated surface as a binary STL file."""
    filename.parent.mkdir(parents=True, exist_ok=True)
    writer = vtk.vtkSTLWriter()
    writer.SetFileName(str(filename))
    writer.SetInputData(surface)
    writer.SetFileTypeToBinary()
    if writer.Write() != 1 or not filename.is_file():
        raise RuntimeError(f"Could not write {filename}")


def compute_centerlines(
    surface: vtk.vtkPolyData, profiles: Sequence[tuple[float, Point]]
) -> vtk.vtkPolyData:
    try:
        from vmtk import vtkvmtk
    except (ImportError, OSError) as exc:
        raise RuntimeError(
            "VMTK and its compatible shared libraries are required for "
            "centerlines. Run this script with the Python executable from the "
            "project Conda environment. Original import error: " + str(exc)
        ) from exc

    # VMTK caps with a new point at each profile centre. These points make
    # stable source/target seeds and avoid interactive endpoint selection.
    capper = vtkvmtk.vtkvmtkCapPolyData()
    capper.SetInputData(surface)
    capper.SetDisplacement(0.0)
    capper.SetInPlaneDisplacement(0.0)
    capper.Update()
    working_surface = capper.GetOutput()
    cap_ids = capper.GetCapCenterIds()
    if cap_ids.GetNumberOfIds() != len(profiles):
        raise RuntimeError("VMTK did not produce one cap centre per boundary.")

    # Match VMTK's cap order to the independently measured boundary profiles.
    available = set(range(cap_ids.GetNumberOfIds()))
    profile_cap_ids: list[int] = []
    for _, centroid in profiles:
        match = min(
            available,
            key=lambda i: vtk.vtkMath.Distance2BetweenPoints(
                centroid, working_surface.GetPoint(cap_ids.GetId(i))
            ),
        )
        profile_cap_ids.append(cap_ids.GetId(match))
        available.remove(match)

    source_index = max(range(len(profiles)), key=lambda i: profiles[i][0])
    source_ids = vtk.vtkIdList()
    source_ids.InsertNextId(profile_cap_ids[source_index])
    target_ids = vtk.vtkIdList()
    for i, point_id in enumerate(profile_cap_ids):
        if i != source_index:
            target_ids.InsertNextId(point_id)

    centerline_filter = vtkvmtk.vtkvmtkPolyDataCenterlines()
    centerline_filter.SetInputData(working_surface)
    centerline_filter.SetSourceSeedIds(source_ids)
    centerline_filter.SetTargetSeedIds(target_ids)
    centerline_filter.SetRadiusArrayName("MaximumInscribedSphereRadius")
    centerline_filter.SetCostFunction("1/R")
    centerline_filter.SetFlipNormals(False)
    centerline_filter.SetAppendEndPointsToCenterlines(True)
    centerline_filter.SetSimplifyVoronoi(False)
    centerline_filter.SetCenterlineResampling(False)
    centerline_filter.Update()
    if centerline_filter.GetOutput().GetNumberOfLines() == 0:
        raise RuntimeError("VMTK produced no centerline paths.")
    output = vtk.vtkPolyData()
    output.DeepCopy(centerline_filter.GetOutput())
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_stl", type=Path, help="open pulmonary-artery STL")
    parser.add_argument(
        "--output-dir", type=Path, default=None, help="default: beside the input STL"
    )
    parser.add_argument(
        "--prefix", default=None, help="output prefix (default: input filename stem)"
    )
    parser.add_argument(
        "--skip-centerlines",
        action="store_true",
        help="write surfaces without requiring VMTK",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_path = args.input_stl.expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    output_dir = (args.output_dir or input_path.parent).expanduser().resolve()
    prefix = args.prefix or input_path.stem

    uncapped = prepare_surface(read_stl(input_path))
    profiles = [profile_area_and_centroid(loop) for loop in boundary_loops(uncapped)]
    if len(profiles) < 2:
        raise RuntimeError(
            f"Expected at least two open boundaries, but found {len(profiles)}."
        )
    capped = cap_surface_vtk(uncapped)

    uncapped_path = output_dir / f"{prefix}_uncapped.vtp"
    capped_path = output_dir / f"{prefix}_capped.stl"
    write_vtp(uncapped, uncapped_path)
    write_stl(capped, capped_path)

    source_index = max(range(len(profiles)), key=lambda i: profiles[i][0])
    print(f"Detected {len(profiles)} open profiles")
    for i, (area, centroid) in enumerate(profiles):
        role = "SOURCE (largest)" if i == source_index else "target"
        print(f"  {i}: area={area:.8g}, centre={centroid}, {role}")
    print(f"Wrote {uncapped_path}")
    print(f"Wrote {capped_path}")

    if not args.skip_centerlines:
        centerlines = compute_centerlines(uncapped, profiles)
        centerlines_path = output_dir / f"{prefix}_centerlines.vtp"
        write_vtp(centerlines, centerlines_path)
        print(f"Wrote {centerlines_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
