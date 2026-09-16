#!/usr/bin/env python3
"""Map Gaussian-smoothed centerline diameters and wall thickness to a surface."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import vtk
from scipy.spatial import cKDTree
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy


def nonnegative_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0.0:
        raise argparse.ArgumentTypeError("must be finite and greater than or equal to 0")
    return number


def percentage(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or not 0.0 < number <= 100.0:
        raise argparse.ArgumentTypeError("must be a finite percentage in (0, 100]")
    return number


def read_surface(filename: Path) -> vtk.vtkPolyData:
    if filename.suffix.lower() == ".stl":
        reader = vtk.vtkSTLReader()
    elif filename.suffix.lower() == ".vtp":
        reader = vtk.vtkXMLPolyDataReader()
    else:
        raise ValueError("input surface must have an .stl or .vtp extension")
    reader.SetFileName(str(filename))
    reader.Update()
    if reader.GetOutput().GetNumberOfPoints() == 0:
        raise RuntimeError(f"No surface points were read from {filename}")
    surface = vtk.vtkPolyData()
    surface.DeepCopy(reader.GetOutput())
    return surface


def read_centerline_csv(
    filename: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    points: list[list[float]] = []
    diameters: list[float] = []
    branches: list[int] = []
    distances: list[float] = []
    with filename.open(newline="", encoding="utf-8") as csv_file:
        reader = csv.DictReader(csv_file)
        required = {"branch", "x", "y", "z", "distance", "diameter"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise RuntimeError(f"CSV is missing columns: {', '.join(sorted(missing))}")
        for row_number, row in enumerate(reader, start=2):
            try:
                point = [float(row[axis]) for axis in ("x", "y", "z")]
                diameter = float(row["diameter"])
                branch = int(row["branch"])
                distance = float(row["distance"])
            except (TypeError, ValueError) as exc:
                raise RuntimeError(f"Invalid numeric data on CSV row {row_number}") from exc
            if not all(math.isfinite(value) for value in (*point, diameter, distance)):
                raise RuntimeError(f"Non-finite numeric data on CSV row {row_number}")
            if diameter <= 0.0:
                raise RuntimeError(f"Non-positive diameter on CSV row {row_number}")
            points.append(point)
            diameters.append(diameter)
            branches.append(branch)
            distances.append(distance)
    if not points:
        raise RuntimeError(f"No valid centerline points were found in {filename}")
    return (
        np.asarray(points, dtype=float),
        np.asarray(diameters, dtype=float),
        np.asarray(branches, dtype=int),
        np.asarray(distances, dtype=float),
    )


def smooth_diameters(
    diameters: np.ndarray,
    branches: np.ndarray,
    distances: np.ndarray,
    sigma: float,
) -> np.ndarray:
    if sigma == 0.0:
        return diameters.copy()
    smoothed = np.empty_like(diameters)
    for branch in np.unique(branches):
        ids = np.flatnonzero(branches == branch)
        branch_distance = distances[ids]
        branch_diameter = diameters[ids]
        # Evaluate the Gaussian in physical distance along the branch, rather
        # than in point indices, so nonuniform centerline sampling is supported.
        for local_index, global_index in enumerate(ids):
            delta = branch_distance - branch_distance[local_index]
            weights = np.exp(-0.5 * (delta / sigma) ** 2)
            weight_sum = weights.sum()
            if weight_sum <= 0.0 or not math.isfinite(float(weight_sum)):
                raise RuntimeError(f"Gaussian weights failed on branch {branch}")
            smoothed[global_index] = np.dot(weights, branch_diameter) / weight_sum
    return smoothed


def nonnegative_int(value: str) -> int:
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError("must be greater than or equal to 0")
    return number


def build_point_adjacency(surface: vtk.vtkPolyData):
    """Unit-weight point adjacency of the surface, its degree, and its edges.

    Coincident points are linked as well, so a surface read from an STL, whose
    triangles do not share point ids, still smooths across triangle boundaries.
    The surface's own points are never modified: mapExtrudeDistance requires the
    written POINT_DATA to have exactly one value per point of this surface.
    """
    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputData(surface)
    triangles.PassLinesOff()
    triangles.PassVertsOff()
    triangles.Update()
    mesh = triangles.GetOutput()
    if mesh.GetNumberOfPoints() != surface.GetNumberOfPoints():
        raise RuntimeError("triangulation changed the surface point count")

    connectivity = vtk_to_numpy(mesh.GetPolys().GetData()).reshape(-1, 4)[:, 1:]
    edges = np.vstack(
        [
            connectivity[:, [0, 1]],
            connectivity[:, [1, 2]],
            connectivity[:, [2, 0]],
        ]
    )

    points = vtk_to_numpy(surface.GetPoints().GetData())
    _, inverse = np.unique(np.round(points, 9), axis=0, return_inverse=True)
    order = np.argsort(inverse, kind="stable")
    same = inverse[order[:-1]] == inverse[order[1:]]
    if same.any():
        edges = np.vstack([edges, np.column_stack([order[:-1][same], order[1:][same]])])

    count = surface.GetNumberOfPoints()
    adjacency = sp.coo_matrix(
        (np.ones(len(edges)), (edges[:, 0], edges[:, 1])), shape=(count, count)
    )
    adjacency = (adjacency + adjacency.T).tocsr()
    adjacency.data[:] = 1.0
    degree = np.maximum(np.asarray(adjacency.sum(1)).ravel(), 1.0)
    return adjacency, degree, edges


def smooth_over_surface(adjacency, degree, values, passes, relaxation=0.5):
    """Laplacian-smooth a scalar field over the surface point graph."""
    for _ in range(passes):
        values = (1.0 - relaxation) * values + relaxation * (adjacency @ values) / degree
    return values


def maximum_field_gradient(points, edges, values) -> float:
    """Largest change in the field per unit distance across any surface edge."""
    length = np.linalg.norm(points[edges[:, 0]] - points[edges[:, 1]], axis=1)
    usable = length > 1e-9
    if not usable.any():
        return 0.0
    difference = np.abs(values[edges[:, 0]] - values[edges[:, 1]])[usable]
    return float((difference / length[usable]).max())


def add_surface_arrays(
    surface: vtk.vtkPolyData,
    centerline_points: np.ndarray,
    smoothed_diameters: np.ndarray,
    extrusion_percentage: float,
    smoothing_passes: int = 10,
) -> tuple[float, float, float, float, float]:
    surface_points = vtk_to_numpy(surface.GetPoints().GetData())
    surface_distances, nearest = cKDTree(centerline_points).query(surface_points)
    mapped_diameters = smoothed_diameters[nearest]

    # The nearest-centerline lookup is discontinuous: two adjacent surface
    # points either side of a bifurcation can snap to different branches, which
    # steps the mapped diameter and so the extrusion distance. varExtrudeMesh
    # offsets each point along its own normal, so a step in the distance twists
    # the extruded cells, and where it coincides with a sharp convex ridge on
    # the wall it produces warped and incorrectly oriented faces in the solid.
    #
    # Smoothing over the surface removes those steps. It is applied to the
    # diameter so that thickness stays exactly the requested percentage of the
    # local diameter, and it flattens spikes without eroding the wall: the
    # minimum, maximum and mean thickness are essentially unchanged.
    adjacency, degree, edges = build_point_adjacency(surface)
    gradient_before = maximum_field_gradient(
        surface_points, edges, mapped_diameters * (extrusion_percentage / 100.0)
    )
    if smoothing_passes > 0:
        mapped_diameters = smooth_over_surface(
            adjacency, degree, mapped_diameters, smoothing_passes
        )

    mapped_thickness = mapped_diameters * (extrusion_percentage / 100.0)
    gradient_after = maximum_field_gradient(surface_points, edges, mapped_thickness)
    diameter_array = numpy_to_vtk(mapped_diameters, deep=True)
    diameter_array.SetName("MappedDiameter")
    thickness_array = numpy_to_vtk(mapped_thickness, deep=True)
    thickness_array.SetName("WallThickness")
    surface.GetPointData().AddArray(diameter_array)
    surface.GetPointData().AddArray(thickness_array)
    surface.GetPointData().SetActiveScalars("WallThickness")
    return (
        float(surface_distances.max()),
        float(mapped_thickness.min()),
        float(mapped_thickness.max()),
        gradient_before,
        gradient_after,
    )


def write_surface(surface: vtk.vtkPolyData, filename: Path) -> None:
    if filename.suffix.lower() != ".vtk":
        raise ValueError("output surface must have a .vtk extension")
    filename.parent.mkdir(parents=True, exist_ok=True)
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(str(filename))
    writer.SetInputData(surface)
    # mapExtrudeDistance reads legacy ASCII VTK POINT_DATA.
    writer.SetFileTypeToASCII()
    if writer.Write() != 1 or not filename.is_file():
        raise RuntimeError(f"Could not write {filename}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_surface", type=Path, help="inner surface STL or VTP")
    parser.add_argument("centerline_csv", type=Path, help="centerline CSV")
    parser.add_argument("output_vtk", type=Path, help="output legacy VTK thickness map")
    parser.add_argument(
        "--gaussian-sigma",
        type=nonnegative_float,
        default=10.0,
        help="diameter smoothing sigma in model units; 0 disables smoothing",
    )
    parser.add_argument(
        "--extrusion-percentage",
        type=percentage,
        default=5.0,
        help="wall thickness as a percentage of local smoothed diameter",
    )
    parser.add_argument(
        "--thickness-smoothing-passes",
        type=nonnegative_int,
        default=10,
        help=(
            "Laplacian passes used to smooth the mapped thickness over the "
            "surface, removing the step the nearest-centerline lookup leaves at "
            "bifurcations; 0 disables"
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_surface = args.input_surface.expanduser().resolve()
    centerline_csv = args.centerline_csv.expanduser().resolve()
    output_vtk = args.output_vtk.expanduser().resolve()
    for filename in (input_surface, centerline_csv):
        if not filename.is_file():
            raise FileNotFoundError(filename)
    surface = read_surface(input_surface)
    points, diameters, branches, distances = read_centerline_csv(centerline_csv)
    smoothed = smooth_diameters(diameters, branches, distances, args.gaussian_sigma)
    (
        maximum_distance,
        minimum_thickness,
        maximum_thickness,
        gradient_before,
        gradient_after,
    ) = add_surface_arrays(
        surface,
        points,
        smoothed,
        args.extrusion_percentage,
        args.thickness_smoothing_passes,
    )
    write_surface(surface, output_vtk)
    print(f"Gaussian sigma: {args.gaussian_sigma:g}")
    print(f"Extrusion percentage: {args.extrusion_percentage:g}%")
    if args.thickness_smoothing_passes > 0:
        print(
            f"Smoothed thickness over the surface: "
            f"{args.thickness_smoothing_passes} passes, peak gradient "
            f"{gradient_before:.4g} to {gradient_after:.4g}"
        )
    print(f"Wall-thickness range: {minimum_thickness:.8g} to {maximum_thickness:.8g}")
    print(f"Maximum surface-to-centerline distance: {maximum_distance:.8g}")
    print(f"Wrote {output_vtk}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, OSError, RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
