#!/usr/bin/env python3
"""Convert an open pulmonary-artery STL to VTP and extract centerlines.

The largest open boundary (by planar profile area) is used as the source;
every other open boundary is used as a target.  VTK performs all surface I/O
and preparation.  VMTK's Python bindings perform the centerline extraction.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import defaultdict, deque
from pathlib import Path
from typing import Sequence

import numpy as np
import scipy.sparse as sp
import vtk
from scipy.spatial import cKDTree
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy


Point = tuple[float, float, float]


def read_stl(filename: Path) -> vtk.vtkPolyData:
    reader = vtk.vtkSTLReader()
    reader.SetFileName(str(filename))
    reader.Update()
    if reader.GetOutput().GetNumberOfCells() == 0:
        raise RuntimeError(f"No triangles could be read from {filename}")
    return reader.GetOutput()


def fill_small_holes(
    surface: vtk.vtkPolyData, minimum_radius: float
) -> tuple[vtk.vtkPolyData, int]:
    """Close open boundaries too small to be a vessel outlet.

    Every open boundary is otherwise promoted to an inlet or outlet: it is
    capped, patched, flow-extended and refined. A pinhole left by a scan defect
    or by a repair tool is therefore treated as a vessel, which both miscounts
    the outlets and asks the mesher for cells orders of magnitude below the
    global size. Holes below minimum_radius are filled into the wall instead.
    """
    if minimum_radius <= 0.0:
        return surface, 0

    before = len(boundary_loops(surface))
    fill = vtk.vtkFillHolesFilter()
    fill.SetInputData(surface)
    fill.SetHoleSize(minimum_radius)
    fill.Update()

    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputConnection(fill.GetOutputPort())
    triangles.PassLinesOff()
    triangles.PassVertsOff()
    triangles.Update()

    # This now runs after prepare_surface's own normals pass, so the lid's
    # winding has to be made consistent here. Splitting stays off to preserve
    # the point count.
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(triangles.GetOutputPort())
    normals.ConsistencyOn()
    normals.AutoOrientNormalsOn()
    normals.SplittingOff()
    normals.Update()

    output = vtk.vtkPolyData()
    output.DeepCopy(normals.GetOutput())
    return output, before - len(boundary_loops(output))


def prepare_surface(
    surface: vtk.vtkPolyData,
    smoothing_iterations: int = 20,
    pass_band: float = 0.05,
) -> vtk.vtkPolyData:
    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(surface)
    clean.PointMergingOn()
    clean.Update()

    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputConnection(clean.GetOutputPort())
    triangles.PassLinesOff()
    triangles.PassVertsOff()
    triangles.Update()

    # Segmentation-derived surfaces carry staircase artefacts whose curvature is
    # far higher than anything anatomical. A point-normal offset folds wherever
    # the wall thickness exceeds the local concave radius of curvature, so this
    # noise seeds self-intersections in the extruded solid. Windowed-sinc
    # (Taubin) smoothing removes it without the shrinkage of a plain Laplacian.
    #
    # Boundary smoothing must stay off: moving the open rims would warp the
    # inlet and outlet profiles that later become fluid patches and the solid's
    # planar symmetry rings. Smoothing precedes the normals filter so the
    # normals describe the final geometry.
    upstream = triangles
    if smoothing_iterations > 0:
        smoother = vtk.vtkWindowedSincPolyDataFilter()
        smoother.SetInputConnection(triangles.GetOutputPort())
        smoother.SetNumberOfIterations(smoothing_iterations)
        smoother.SetPassBand(pass_band)
        smoother.BoundarySmoothingOff()
        smoother.FeatureEdgeSmoothingOff()
        smoother.NonManifoldSmoothingOn()
        smoother.NormalizeCoordinatesOn()
        smoother.Update()
        upstream = smoother

    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(upstream.GetOutputPort())
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


def _point_adjacency(surface: vtk.vtkPolyData):
    """Unit-weight point adjacency of a triangulated surface, and its degree."""
    polygons = vtk_to_numpy(surface.GetPolys().GetData()).reshape(-1, 4)[:, 1:]
    edges = np.vstack(
        [polygons[:, [0, 1]], polygons[:, [1, 2]], polygons[:, [2, 0]]]
    )
    count = surface.GetNumberOfPoints()
    adjacency = sp.coo_matrix(
        (np.ones(len(edges)), (edges[:, 0], edges[:, 1])), shape=(count, count)
    )
    adjacency = (adjacency + adjacency.T).tocsr()
    adjacency.data[:] = 1.0
    degree = np.maximum(np.asarray(adjacency.sum(1)).ravel(), 1.0)
    return adjacency, degree


def _smooth_scalar(adjacency, degree, values, passes, relaxation=0.5):
    for _ in range(passes):
        values = (1.0 - relaxation) * values + relaxation * (adjacency @ values) / degree
    return values


def _dilate_scalar(adjacency, values, passes):
    """Grow a mask by repeated neighbourhood maximum."""
    for _ in range(passes):
        values = np.maximum(
            values,
            np.asarray(
                sp.csr_matrix(adjacency.multiply(values[None, :])).max(axis=1).todense()
            ).ravel(),
        )
    return values


def _minimum_principal_curvature(surface: vtk.vtkPolyData) -> np.ndarray:
    """Most concave principal curvature at every point."""
    gaussian = vtk.vtkCurvatures()
    gaussian.SetInputData(surface)
    gaussian.SetCurvatureTypeToGaussian()
    gaussian.Update()
    mean = vtk.vtkCurvatures()
    mean.SetInputData(surface)
    mean.SetCurvatureTypeToMean()
    mean.Update()
    K = vtk_to_numpy(gaussian.GetOutput().GetPointData().GetArray("Gauss_Curvature"))
    H = vtk_to_numpy(mean.GetOutput().GetPointData().GetArray("Mean_Curvature"))
    return H - np.sqrt(np.maximum(H * H - K, 0.0))


def centerline_branch_arrays(centerlines: vtk.vtkPolyData):
    """Per-point centerline coordinates, radii, branch index and arc length.

    This mirrors the traversal in centerlines_to_csv.py so the diameters used
    here are the same ones Gauss.py will later turn into wall thickness.
    """
    radius_array = centerlines.GetPointData().GetArray("MaximumInscribedSphereRadius")
    if radius_array is None:
        raise RuntimeError(
            "Centerlines are missing the MaximumInscribedSphereRadius array"
        )

    coordinates: list[Point] = []
    radii: list[float] = []
    branches: list[int] = []
    distances: list[float] = []

    lines = centerlines.GetLines()
    lines.InitTraversal()
    point_ids = vtk.vtkIdList()
    branch = 0
    while lines.GetNextCell(point_ids):
        distance = 0.0
        previous: Point | None = None
        for index in range(point_ids.GetNumberOfIds()):
            point_id = point_ids.GetId(index)
            point = centerlines.GetPoint(point_id)
            if previous is not None:
                distance += math.dist(previous, point)
            coordinates.append(point)
            radii.append(radius_array.GetTuple1(point_id))
            branches.append(branch)
            distances.append(distance)
            previous = point
        branch += 1

    return (
        np.asarray(coordinates, dtype=float),
        np.asarray(radii, dtype=float),
        np.asarray(branches, dtype=int),
        np.asarray(distances, dtype=float),
    )


def local_wall_thickness(
    surface: vtk.vtkPolyData,
    centerlines: vtk.vtkPolyData,
    extrusion_percentage: float,
    gaussian_sigma: float,
) -> np.ndarray:
    """Wall thickness that will later be extruded at each surface point.

    Gauss.py's own smoothing and nearest-centerline mapping are reused so the
    threshold here matches the thickness that is actually extruded.
    """
    from Gauss import smooth_diameters

    points, radii, branches, distances = centerline_branch_arrays(centerlines)
    smoothed = smooth_diameters(2.0 * radii, branches, distances, gaussian_sigma)
    surface_points = vtk_to_numpy(surface.GetPoints().GetData())
    _, nearest = cKDTree(points).query(surface_points)
    return smoothed[nearest] * (extrusion_percentage / 100.0)


def _surface_triangles(surface: vtk.vtkPolyData) -> np.ndarray:
    """Triangle connectivity of the surface, in its own point numbering."""
    filter_ = vtk.vtkTriangleFilter()
    filter_.SetInputData(surface)
    filter_.PassLinesOff()
    filter_.PassVertsOff()
    filter_.Update()
    mesh = filter_.GetOutput()
    if mesh.GetNumberOfPoints() != surface.GetNumberOfPoints():
        raise RuntimeError("triangulation changed the surface point count")
    return vtk_to_numpy(mesh.GetPolys().GetData()).reshape(-1, 4)[:, 1:]


def _triangle_areas(triangles: np.ndarray, points: np.ndarray) -> np.ndarray:
    a, b, c = points[triangles[:, 0]], points[triangles[:, 1]], points[triangles[:, 2]]
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)


def _collapsing_triangles(
    triangles: np.ndarray, initial_areas: np.ndarray, candidate: np.ndarray
) -> np.ndarray:
    """Triangles that candidate would shrink far below their starting area.

    The comparison is against the area before any smoothing, not the previous
    pass: a per-pass limit still permits unbounded collapse once it is applied a
    few hundred times.
    """
    return _triangle_areas(triangles, candidate) < 0.25 * np.maximum(
        initial_areas, 1e-30
    )


def smooth_saddles(
    surface: vtk.vtkPolyData,
    target_radius,
    dilation: int = 3,
    outer_iterations: int = 20,
) -> tuple[vtk.vtkPolyData, float]:
    """Round only the concave saddles until their radius of curvature grows.

    A point-normal offset folds where the wall thickness exceeds the local
    concave radius of curvature, and on a vessel tree that happens almost
    exclusively in the bifurcation crotches. Global smoothing cannot fix this
    without deforming the whole anatomy, so this smooths only the points whose
    concave radius is below target_radius, with the region grown and feathered
    so the correction blends into untouched surface. It is the automated form
    of filleting the crotches by hand.

    target_radius is either one value for the whole surface or an array with
    one value per surface point, which lets the threshold follow the wall
    thickness that will actually be extruded there.

    Open-profile rim points are held fixed. Returns the smoothed surface and
    the fraction of points that were moved appreciably.
    """
    working = vtk.vtkPolyData()
    working.DeepCopy(surface)
    adjacency, degree = _point_adjacency(working)
    triangles = _surface_triangles(working)
    original = vtk_to_numpy(working.GetPoints().GetData()).copy()
    initial_areas = _triangle_areas(triangles, original)

    # The rims bound the fluid patches and the solid's planar symmetry rings.
    pinned = np.zeros(working.GetNumberOfPoints(), dtype=bool)
    edges = vtk.vtkFeatureEdges()
    edges.SetInputData(working)
    edges.BoundaryEdgesOn()
    edges.FeatureEdgesOff()
    edges.NonManifoldEdgesOff()
    edges.ManifoldEdgesOff()
    edges.Update()
    locator = vtk.vtkPointLocator()
    locator.SetDataSet(working)
    locator.BuildLocator()
    rim = edges.GetOutput()
    for i in range(rim.GetNumberOfPoints()):
        pinned[locator.FindClosestPoint(rim.GetPoint(i))] = True

    for _ in range(outer_iterations):
        curvature = _smooth_scalar(
            adjacency, degree, _minimum_principal_curvature(working), 3
        )
        radius = np.where(
            curvature < 0.0, 1.0 / np.maximum(np.abs(curvature), 1e-9), np.inf
        )
        target = np.maximum(np.asarray(target_radius, dtype=float), 1e-9)
        weight = np.clip((target - radius) / target, 0.0, 1.0)
        weight = _dilate_scalar(adjacency, weight, dilation)
        weight = _smooth_scalar(adjacency, degree, weight, 5)
        weight[pinned] = 0.0

        points = vtk_to_numpy(working.GetPoints().GetData())
        for _ in range(10):
            candidate = points + 0.6 * weight[:, None] * (
                (adjacency @ points) / degree[:, None] - points
            )

            # Hold back any point whose movement would collapse one of its
            # triangles. Smoothing can drive three points collinear while every
            # edge stays long, and a zero-area triangle has no usable normal
            # exactly where the extrusion needs one; it also makes the thickness
            # map unprojectable onto the wall patch.
            collapsing = _collapsing_triangles(triangles, initial_areas, candidate)
            if collapsing.any():
                frozen = np.zeros(len(points), dtype=bool)
                frozen[triangles[collapsing].ravel()] = True
                candidate[frozen] = points[frozen]
            points = candidate
        working.GetPoints().SetData(numpy_to_vtk(points, deep=True))
        working.Modified()

    moved = np.linalg.norm(
        vtk_to_numpy(working.GetPoints().GetData()) - original, axis=1
    )
    return working, float((moved > 0.02).sum()) / len(moved)


def import_vmtk():
    """Import VMTK's VTK bindings with an actionable error message."""
    try:
        from vmtk import vtkvmtk
    except (ImportError, OSError) as exc:
        raise RuntimeError(
            "VMTK and its compatible shared libraries are required. Run this "
            "script with the Python executable from the project Conda "
            "environment. Original import error: " + str(exc)
        ) from exc
    return vtkvmtk


# Fraction of each extension over which the real, generally non-circular
# profile is blended into the circular rim VMTK always produces. Measured on
# these geometries: 0.5 leaves a 63 degree crease that the solid extrusion
# folds, 0.8 leaves 41 degrees, and going further buys nothing.
EXTENSION_TRANSITION_RATIO = 0.8


def extend_profiles(
    surface: vtk.vtkPolyData,
    centerlines: vtk.vtkPolyData,
    extension_diameters: float,
) -> vtk.vtkPolyData:
    """Add straight flow extensions of the given length to every open profile.

    Each extension follows the local centerline direction and terminates in a
    circular, planar rim. This makes capping exact, moves the inlet and outlet
    boundary conditions off the bifurcations, and gives the extruded solid end
    rings that are genuinely planar and normal to the vessel axis, as the
    symmetry patches assigned to them require.
    """
    vtkvmtk = import_vmtk()

    extensions = vtkvmtk.vtkvmtkPolyDataFlowExtensionsFilter()
    extensions.SetInputData(surface)
    extensions.SetCenterlines(centerlines)

    # Extend along the vessel axis rather than the rim normal: on an obliquely
    # cut profile the rim normal is itself oblique.
    extensions.SetExtensionModeToUseCenterlineDirection()
    extensions.SetInterpolationModeToThinPlateSpline()

    # With adaptive length the extension is ExtensionRatio times the profile's
    # mean radius, so one diameter corresponds to a ratio of two.
    extensions.SetAdaptiveExtensionLength(1)
    extensions.SetAdaptiveExtensionRadius(1)
    extensions.SetAdaptiveNumberOfBoundaryPoints(0)
    extensions.SetExtensionRatio(2.0 * extension_diameters)
    extensions.SetTransitionRatio(EXTENSION_TRANSITION_RATIO)
    extensions.SetCenterlineNormalEstimationDistanceRatio(1.0)
    extensions.SetNumberOfBoundaryPoints(50)
    extensions.SetSigma(1.0)
    extensions.Update()

    # The extension filter emits triangles whose winding does not follow the
    # surface it was given, which leaves pairs of neighbouring triangles with
    # opposing normals. Nothing downstream re-orients them: prepare_surface runs
    # before this stage. Left alone they reach cfMesh as apparent creases and
    # end up as warped, incorrectly oriented faces in the extruded solid.
    # Splitting must stay off so the point count is preserved.
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(extensions.GetOutputPort())
    normals.ConsistencyOn()
    normals.AutoOrientNormalsOn()
    normals.SplittingOff()
    normals.Update()

    output = vtk.vtkPolyData()
    output.DeepCopy(normals.GetOutput())
    if output.GetNumberOfCells() == 0:
        raise RuntimeError("Flow extension produced an empty surface.")
    return output


CELL_ENTITY_IDS = "CellEntityIds"


def cap_surface(surface: vtk.vtkPolyData) -> tuple[vtk.vtkPolyData, np.ndarray]:
    """Close every open profile, tagging the wall and each cap separately.

    The capper labels the cells it adds, which is what lets every opening be
    named here rather than recovered later from the volume mesh. Deducing the
    openings after meshing meant separating boundary faces by the angle between
    them, and that cannot tell a cap from a patch of wall: a region refined
    finer than its surroundings was split off as though it were an opening,
    while a cap only a few cells across was absorbed into the wall.

    Returns the capped surface and one entity id per cell.
    """
    vtkvmtk = import_vmtk()

    capper = vtkvmtk.vtkvmtkCapPolyData()
    capper.SetInputData(surface)
    capper.SetDisplacement(0.0)
    capper.SetInPlaneDisplacement(0.0)
    capper.SetCellEntityIdsArrayName(CELL_ENTITY_IDS)
    capper.SetCellEntityIdOffset(0)
    capper.Update()

    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputConnection(capper.GetOutputPort())
    triangles.Update()
    output = vtk.vtkPolyData()
    output.DeepCopy(triangles.GetOutput())

    array = output.GetCellData().GetArray(CELL_ENTITY_IDS)
    if array is None:
        raise RuntimeError("the capper did not tag its cells with entity ids")
    return output, vtk_to_numpy(array).astype(int)


def name_cap_regions(
    surface: vtk.vtkPolyData, entity_ids: np.ndarray, profiles
) -> dict[int, str]:
    """Give each entity id the patch name it should carry into the mesh.

    The wall is the region with by far the most cells. Each remaining region is
    matched to the open profile it sits on, and the one on the largest profile
    becomes the inlet, as that is the main pulmonary artery.
    """
    points = vtk_to_numpy(surface.GetPoints().GetData())
    polygons = vtk_to_numpy(surface.GetPolys().GetData()).reshape(-1, 4)[:, 1:]
    centres = np.array([centroid for _, centroid in profiles])
    areas = np.array([area for area, _ in profiles])

    unique, counts = np.unique(entity_ids, return_counts=True)
    wall_id = int(unique[counts.argmax()])

    matched: list[tuple[int, int]] = []
    for identifier in unique:
        if int(identifier) == wall_id:
            continue
        selection = polygons[entity_ids == identifier]
        centre = points[np.unique(selection)].mean(0)
        matched.append(
            (int(identifier), int(np.linalg.norm(centres - centre, axis=1).argmin()))
        )

    if not matched:
        raise RuntimeError("the capper produced no caps")

    inlet_id = max(matched, key=lambda pair: areas[pair[1]])[0]
    names = {wall_id: "wall", inlet_id: "inlet"}
    for number, (identifier, _) in enumerate(
        (pair for pair in matched if pair[0] != inlet_id), start=1
    ):
        names[identifier] = f"outlet{number}"
    return names


TARGET_AREA = "TargetArea"


def remesh_caps(
    surface: vtk.vtkPolyData, entity_ids: np.ndarray, wall_id: int
) -> tuple[vtk.vtkPolyData, np.ndarray]:
    """Replace each cap's triangle fan with an isotropic triangulation.

    The capper resamples a rim into a regular polygon and fans it to a single
    centre point, so every cap arrives as a ring of congruent slivers: on p02
    each cap was 50 identical triangles with a 7.2 degree apex, a normalised
    shape quality of 0.215 against 1.0 for an equilateral triangle. That is one
    long triangle spanning the whole opening in place of a mesh, which leaves
    cfMesh nothing to project a boundary vertex onto across the cap interior and
    gives surfaceFeatureEdges badly conditioned triangles to take dihedral
    angles from.

    Only the caps are rebuilt. The wall is handed to the remesher as an excluded
    region, so its triangles are not touched and the rim shared with each cap
    cannot move; the vessel surface that the thickness mapping and the solid
    extrusion depend on is therefore bit-for-bit what it was. The caps stay in
    their own planes, which the extensions made planar, so the openings remain
    flat.

    Each cap is sized from the wall it meets rather than from a single global
    length, because the two ends of one geometry are nothing like each other: on
    p02 the wall triangles around the inlet average 1.87 mm and those around
    outlet11 average 0.187 mm. A shared target would leave the smallest outlets
    with triangles wider than the outlet itself.
    """
    vtkvmtk = import_vmtk()

    points = vtk_to_numpy(surface.GetPoints().GetData())
    polygons = vtk_to_numpy(surface.GetPolys().GetData()).reshape(-1, 4)[:, 1:]
    wall_triangles = polygons[entity_ids == wall_id]
    if not len(wall_triangles):
        raise RuntimeError("no wall cells were found to remesh the caps against")

    def mean_edge_length(triangles: np.ndarray) -> float:
        corners = points[triangles]
        return float(np.linalg.norm(corners - np.roll(corners, 1, 1), axis=2).mean())

    # The target is point data, and a rim point belongs to both its cap and the
    # wall. Giving it the cap's target is what is wanted: the new triangles die
    # down to the size of the wall triangles they meet.
    targets = np.full(len(points), mean_edge_length(wall_triangles))
    for identifier in np.unique(entity_ids):
        if identifier == wall_id:
            continue
        cap_points = np.unique(polygons[entity_ids == identifier])
        neighbours = wall_triangles[np.isin(wall_triangles, cap_points).any(1)]
        if not len(neighbours):
            raise RuntimeError(f"cap {identifier} does not touch the wall")
        targets[cap_points] = mean_edge_length(neighbours)

    working = vtk.vtkPolyData()
    working.DeepCopy(surface)
    # The capper writes the entity ids as a vtkIdTypeArray but the remesher
    # reads a vtkIntArray. Left as it is, the array is silently not found, no
    # region is excluded and the whole surface is remeshed.
    identifiers = numpy_to_vtk(
        entity_ids.astype(np.int32), deep=1, array_type=vtk.VTK_INT
    )
    identifiers.SetName(CELL_ENTITY_IDS)
    working.GetCellData().AddArray(identifiers)
    areas = numpy_to_vtk(0.25 * math.sqrt(3.0) * targets**2, deep=1)
    areas.SetName(TARGET_AREA)
    working.GetPointData().AddArray(areas)

    excluded = vtk.vtkIdList()
    excluded.InsertNextId(int(wall_id))

    remesher = vtkvmtk.vtkvmtkPolyDataSurfaceRemeshing()
    remesher.SetInputData(working)
    remesher.SetCellEntityIdsArrayName(CELL_ENTITY_IDS)
    remesher.SetExcludedEntityIds(excluded)
    remesher.SetElementSizeModeToTargetAreaArray()
    remesher.SetTargetAreaArrayName(TARGET_AREA)
    remesher.SetNumberOfIterations(10)
    remesher.SetNumberOfConnectivityOptimizationIterations(20)
    remesher.SetPreserveBoundaryEdges(1)
    # The filter counts its iterations out on standard output from C++, which
    # tee'd to the screen buries the one line this stage is meant to print.
    # Errors go to standard error and are not affected.
    sys.stdout.flush()
    saved = os.dup(1)
    try:
        with open(os.devnull, "w") as quiet:
            os.dup2(quiet.fileno(), 1)
        remesher.Update()
    finally:
        os.dup2(saved, 1)
        os.close(saved)

    output = vtk.vtkPolyData()
    output.DeepCopy(remesher.GetOutput())
    array = output.GetCellData().GetArray(CELL_ENTITY_IDS)
    if array is None:
        raise RuntimeError("the remesher did not carry the entity ids through")
    updated = vtk_to_numpy(array).astype(int)

    # Every patch must survive, and the wall must come back untouched: if the
    # entity ids were not read the wall would have been remeshed too, silently
    # moving the surface the solid is extruded from.
    if set(updated.tolist()) != set(entity_ids.tolist()):
        raise RuntimeError("cap remeshing lost or invented a named region")
    before = np.sort(points[wall_triangles].reshape(-1, 9), axis=0)
    new_points = vtk_to_numpy(output.GetPoints().GetData())
    new_polygons = vtk_to_numpy(output.GetPolys().GetData()).reshape(-1, 4)[:, 1:]
    after = np.sort(new_points[new_polygons[updated == wall_id]].reshape(-1, 9), axis=0)
    if before.shape != after.shape or not np.allclose(before, after):
        raise RuntimeError("cap remeshing moved the vessel wall")

    edges = vtk.vtkFeatureEdges()
    edges.SetInputData(output)
    edges.BoundaryEdgesOn()
    edges.NonManifoldEdgesOn()
    edges.FeatureEdgesOff()
    edges.ManifoldEdgesOff()
    edges.Update()
    if edges.GetOutput().GetNumberOfCells():
        raise RuntimeError("cap remeshing left the surface open or non-manifold")

    return output, updated


def write_named_stl(
    surface: vtk.vtkPolyData,
    entity_ids: np.ndarray,
    names: dict[int, str],
    filename: Path,
) -> None:
    """Write one ASCII STL solid per named region.

    cfMesh turns each solid into a patch of the same name, so the patches of the
    finished mesh are decided here, by construction, rather than inferred.
    """
    filename.parent.mkdir(parents=True, exist_ok=True)
    points = vtk_to_numpy(surface.GetPoints().GetData())
    polygons = vtk_to_numpy(surface.GetPolys().GetData()).reshape(-1, 4)[:, 1:]
    a, b, c = points[polygons[:, 0]], points[polygons[:, 1]], points[polygons[:, 2]]
    normals = np.cross(b - a, c - a)
    normals /= np.maximum(np.linalg.norm(normals, axis=1), 1e-30)[:, None]

    with filename.open("w", encoding="ascii") as handle:
        for identifier, name in sorted(names.items(), key=lambda kv: kv[1]):
            handle.write(f"solid {name}\n")
            for index in np.flatnonzero(entity_ids == identifier):
                handle.write(
                    " facet normal %.9g %.9g %.9g\n  outer loop\n"
                    % tuple(normals[index])
                )
                for vertex in (a[index], b[index], c[index]):
                    handle.write("   vertex %.9g %.9g %.9g\n" % tuple(vertex))
                handle.write("  endloop\n endfacet\n")
            handle.write(f"endsolid {name}\n")
    if not filename.is_file():
        raise RuntimeError(f"Could not write {filename}")


def write_profiles_csv(profiles, filename: Path) -> None:
    """Record each open profile's size and centre for the meshing stage.

    These are the rims of the surface that is actually meshed, so with flow
    extensions enabled they are the extended rims, not the original ones.
    """
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["index", "area", "radius", "x", "y", "z"])
        for index, (area, centroid) in enumerate(profiles):
            radius = math.sqrt(area / math.pi) if area > 0.0 else 0.0
            writer.writerow(
                [index, area, radius, centroid[0], centroid[1], centroid[2]]
            )


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
    vtkvmtk = import_vmtk()

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
        help="skip centerline extraction, and with it the flow extensions",
    )
    parser.add_argument(
        "--minimum-profile-radius",
        type=float,
        default=1.0,
        help=(
            "open boundaries smaller than this radius, in model units, are "
            "treated as surface defects and filled rather than being taken for "
            "vessel outlets; 0 keeps every hole"
        ),
    )
    parser.add_argument(
        "--smoothing-iterations",
        type=int,
        default=20,
        help="windowed-sinc surface smoothing iterations; 0 disables smoothing",
    )
    parser.add_argument(
        "--smoothing-passband",
        type=float,
        default=0.05,
        help="windowed-sinc pass band; lower smooths more (default: 0.05)",
    )
    parser.add_argument(
        "--no-saddle-smoothing",
        action="store_true",
        help="do not round concave saddles at all",
    )
    parser.add_argument(
        "--saddle-safety-factor",
        type=float,
        default=0.6,
        help=(
            "a saddle is rounded until its concave radius reaches the local wall "
            "thickness divided by this factor (default: 0.6)"
        ),
    )
    parser.add_argument(
        "--gaussian-sigma",
        type=float,
        default=10.0,
        help="diameter smoothing sigma, matching Gauss.py, for the saddle target",
    )
    parser.add_argument(
        "--extrusion-percentage",
        type=float,
        default=5.0,
        help="wall thickness as a percentage of diameter, matching Gauss.py",
    )
    parser.add_argument(
        "--saddle-dilation",
        type=int,
        default=3,
        help="rings by which the saddle region is grown before feathering",
    )
    parser.add_argument(
        "--extension-diameters",
        type=float,
        default=1.0,
        help="flow-extension length per profile, in local diameters; 0 disables",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_path = args.input_stl.expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    output_dir = (args.output_dir or input_path.parent).expanduser().resolve()
    prefix = args.prefix or input_path.stem

    uncapped = prepare_surface(
        read_stl(input_path),
        smoothing_iterations=args.smoothing_iterations,
        pass_band=args.smoothing_passband,
    )

    # Fill after smoothing, not before. BoundarySmoothingOff protects the points
    # of an open rim, but once a pinhole is closed its lid is ordinary interior
    # surface, and since the hole is smaller than one triangle the smoothing
    # collapses that lid to zero area.
    uncapped, filled = fill_small_holes(uncapped, args.minimum_profile_radius)
    if filled:
        print(
            f"Filled {filled} open boundary/boundaries smaller than "
            f"{args.minimum_profile_radius:g}: too small to be a vessel outlet"
        )
    if args.smoothing_iterations > 0:
        print(
            f"Smoothed the surface: {args.smoothing_iterations} windowed-sinc "
            f"iterations at pass band {args.smoothing_passband:g}"
        )

    profiles = [profile_area_and_centroid(loop) for loop in boundary_loops(uncapped)]
    if len(profiles) < 2:
        raise RuntimeError(
            f"Expected at least two open boundaries, but found {len(profiles)}."
        )

    source_index = max(range(len(profiles)), key=lambda i: profiles[i][0])
    print(f"Detected {len(profiles)} open profiles")
    for i, (area, centroid) in enumerate(profiles):
        role = "SOURCE (largest)" if i == source_index else "target"
        print(f"  {i}: area={area:.8g}, centre={centroid}, {role}")

    # Centerlines describe the original anatomy, so they are extracted before
    # the profiles are extended. They also supply the extension directions.
    # Surface points on an extension map to the nearest centerline endpoint,
    # which carries a constant diameter along the straight extension.
    centerlines = None
    centerlines_path = output_dir / f"{prefix}_centerlines.vtp"
    if not args.skip_centerlines:
        centerlines = compute_centerlines(uncapped, profiles)
        write_vtp(centerlines, centerlines_path)
        print(f"Wrote {centerlines_path}")

    # Round the concave saddles now that the centerlines are available: the
    # threshold at each point is the wall thickness that will be extruded
    # there, divided by the safety factor. Because this moves the surface, the
    # centerlines are re-extracted afterwards so every downstream diameter
    # describes the geometry that is actually meshed.
    if centerlines is not None and not args.no_saddle_smoothing:
        thickness = local_wall_thickness(
            uncapped, centerlines, args.extrusion_percentage, args.gaussian_sigma
        )
        target = thickness / args.saddle_safety_factor
        description = (
            f"the local wall thickness / {args.saddle_safety_factor:g} "
            f"(target radius {target.min():.3g} to {target.max():.3g})"
        )

        uncapped, touched = smooth_saddles(uncapped, target, args.saddle_dilation)
        print(
            f"Rounded concave saddles to {description}: "
            f"{100.0 * touched:.2f}% of the surface was moved"
        )

        if touched > 0.0:
            centerlines = compute_centerlines(uncapped, profiles)
            write_vtp(centerlines, centerlines_path)
            print(f"Re-extracted centerlines on the rounded surface")
    elif not args.no_saddle_smoothing:
        print("Skipping saddle smoothing because centerlines were not extracted.")

    if centerlines is not None and args.extension_diameters > 0.0:
        uncapped = extend_profiles(uncapped, centerlines, args.extension_diameters)
        print(
            f"Extended every profile by {args.extension_diameters:g} diameter(s)"
        )
    elif args.extension_diameters > 0.0:
        print("Skipping flow extensions because centerlines were not extracted.")

    # Re-measure the rims: with flow extensions these are the extended profiles,
    # which is where the caps and therefore the outlet patches end up.
    final_profiles = [
        profile_area_and_centroid(loop) for loop in boundary_loops(uncapped)
    ]
    profiles_path = output_dir / f"{prefix}_profiles.csv"
    write_profiles_csv(final_profiles, profiles_path)
    print(f"Wrote {profiles_path}")

    capped, entity_ids = cap_surface(uncapped)
    names = name_cap_regions(capped, entity_ids, final_profiles)

    wall_id = next(key for key, name in names.items() if name == "wall")
    cap_cells = int((entity_ids != wall_id).sum())
    capped, entity_ids = remesh_caps(capped, entity_ids, wall_id)
    print(
        f"Remeshed the caps: {cap_cells} triangles to "
        f"{int((entity_ids != wall_id).sum())}, the wall untouched"
    )

    uncapped_path = output_dir / f"{prefix}_uncapped.vtp"
    capped_path = output_dir / f"{prefix}_capped.stl"
    write_vtp(uncapped, uncapped_path)
    write_named_stl(capped, entity_ids, names, capped_path)
    outlets = sum(1 for name in names.values() if name.startswith("outlet"))
    print(f"Named the surface: wall, inlet and {outlets} outlet(s)")
    print(f"Wrote {uncapped_path}")
    print(f"Wrote {capped_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
