# mapExtrudeDistance

`mapExtrudeDistance` is an OpenFOAM utility that maps a scalar field from the
`POINT_DATA` section of a legacy ASCII VTK surface onto the points of an
OpenFOAM boundary patch. It writes the result as a `pointScalarField` in the
selected case's current time directory.

Despite the name, this utility **does not extrude the mesh or move any mesh
points**. It only creates a mapped point field that can be used by a later
workflow.

## How the mapping works

For each point on the target boundary patch, the utility:

1. Finds the nearest triangle on the VTK surface using an indexed octree.
2. Projects the target point onto that triangle.
3. Barycentrically interpolates the scalar values stored at the triangle's
   three vertices.
4. Rejects the mapping if the closest source position is farther away than
   `maxDistance`.

Polygonal VTK faces are triangulated internally using a triangle fan. The
source VTK geometry and the OpenFOAM patch must therefore use the same
coordinate system and units.

## Requirements

- A configured OpenFOAM development environment with `wmake` available.
- An existing OpenFOAM case containing a volume mesh.
- A legacy ASCII VTK surface containing polygon geometry and a scalar
  `POINT_DATA` array.

The small field parser does not support binary legacy VTK, XML VTK/VTU,
vector arrays, or `FIELD` arrays.

## Build

Load the appropriate OpenFOAM environment, then run:

```bash
cd mapExtrudeDistance
wmake
```

The executable is installed as:

```text
$FOAM_USER_APPBIN/mapExtrudeDistance
```

The build links against OpenFOAM's `finiteVolume`, `surfMesh`, and `meshTools`
libraries.

## Input VTK format

The input must be a legacy ASCII VTK file. The requested scalar array must
appear in `POINT_DATA` in the standard `SCALARS` form, for example:

```text
# vtk DataFile Version 3.0
Wall thickness
ASCII
DATASET POLYDATA
POINTS 3 float
0 0 0
1 0 0
0 1 0
POLYGONS 1 4
3 0 1 2
POINT_DATA 3
SCALARS WallThickness float 1
LOOKUP_TABLE default
0.001
0.0015
0.002
```

There must be exactly one finite scalar value per imported VTK surface point.
The parser can skip preceding scalar arrays and select the one named by
`-vtkField`.

## Usage

From an OpenFOAM case directory:

```bash
mapExtrudeDistance
```

With explicit settings:

```bash
mapExtrudeDistance \
    -case /path/to/case \
    -vtk constant/triSurface/thickMap.vtk \
    -patch wall \
    -vtkField WallThickness \
    -outputField WallThickness \
    -maxDistance 0.2
```

Relative VTK paths are resolved relative to the selected case directory,
including when `-case` is used.

### Options

| Option | Default | Description |
| --- | --- | --- |
| `-vtk <file>` | `thickMap.vtk` in the case directory | Source legacy ASCII VTK surface. |
| `-patch <name>` | `wall` | Target OpenFOAM boundary patch. |
| `-vtkField <name>` | `WallThickness` | Scalar array to read from VTK `POINT_DATA`. |
| `-outputField <name>` | Same as `-vtkField` | Name of the written OpenFOAM point field. |
| `-maxDistance <scalar>` | `0.2` | Maximum permitted distance from a target point to the source surface. Must be positive and use the geometry's length units. |

Standard OpenFOAM command-line options, such as `-case`, are also available.
Use `mapExtrudeDistance -help` to see the options supported by the installed
build.

## Output

The utility writes a length-dimensioned `pointScalarField` to:

```text
<case>/<current-time>/<outputField>
```

Values are assigned to global mesh points belonging to the selected patch.
All other mesh-point values remain zero. If a mesh point is shared with
another patch, its global primitive value receives the mapped value and the
point-patch fields are then updated by boundary-condition correction.

On completion, the utility reports:

- the number of mapped target points;
- the largest source-surface projection distance;
- the minimum, maximum, and average mapped values; and
- the path of the written field.

## Validation and common errors

The program stops with a fatal error when:

- the requested boundary patch does not exist;
- the VTK file cannot be opened or is not legacy ASCII VTK;
- no `POINT_DATA` section or requested scalar array is found;
- the VTK point-data count differs from the imported surface point count;
- the selected array is not scalar or contains missing/non-finite values;
- the source surface is empty or contains invalid faces;
- `maxDistance` is not positive; or
- any target point has no source triangle within `maxDistance`.

If the distance check fails, first verify that both geometries use the same
coordinates and units. Increase `-maxDistance` only when the separation is
expected and physically appropriate; the check helps prevent silently mapping
from an unrelated part of the surface.

## Source layout

- `mapExtrudeDistance.C` provides the command-line utility and writes the
  OpenFOAM point field.
- `variableDistance.H` and `variableDistance.C` provide the VTK reader and
  surface-to-point interpolation routines.
- `Make/files` and `Make/options` define the OpenFOAM `wmake` build.

The `variableDistance` helper also contains constructors for converting a
strictly positive, face-based distance list into area-weighted point values.
That API is available to other C++ code but is not invoked by the
`mapExtrudeDistance` executable described above.
