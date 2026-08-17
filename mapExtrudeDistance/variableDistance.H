/*---------------------------------------------------------------------------*\
  variableDistance.H

  Purpose
      Declare a small helper class for variable-distance surface extrusion.

      The input distance is defined once per surface face (VTU CellData is a
      natural source for this).  OpenFOAM extrusion ultimately moves shared
      surface points, so the class converts those face values into one value
      per surface point using area-weighted averaging.

  Important
      This header declares the public interface and the stored data.  The
      corresponding function bodies are in variableDistance.C.
\*---------------------------------------------------------------------------*/

#ifndef variableDistance_H
#define variableDistance_H

// OpenFOAM field of scalar values (normally double precision).
#include "scalarField.H"

// OpenFOAM polygonal surface class used by constructFrom surface.
#include "MeshedSurface.H"

// OpenFOAM's general polygon face type.
#include "face.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class variableDistance
\*---------------------------------------------------------------------------*/

class variableDistance
{
private:

    // One extrusion distance for every input surface face.
    //
    // Entry i must refer to surface face i.  This ordering requirement is
    // why geometry and VTU CellData must be exported without independently
    // reordering, cleaning, or triangulating one of them.
    scalarField faceDistance_;


    // One extrusion distance for every input surface point.
    //
    // Adjacent faces share points and cannot move a shared point by different
    // distances.  We therefore calculate this field from faceDistance_ using
    // an area-weighted average over all faces incident on each point.
    scalarField pointDistance_;


    // Private Member Functions

        // Read a raw OpenFOAM scalar list from a file.
        //
        // The expected file intentionally has no FoamFile header.  Example:
        //
        //     4
        //     (
        //         0.001
        //         0.0015
        //         0.002
        //         0.0025
        //     )
        //
        // Keeping this reader simple makes the face ordering completely
        // explicit.  A later VTU reader can supply a scalarField directly by
        // using the second constructor below.
        static scalarField readFaceDistance(const fileName& distanceFile);


        // Read a named scalar array from POINT_DATA in a legacy ASCII VTK
        // file.  The parser supports the common form
        //
        //     POINT_DATA <numberOfPoints>
        //     SCALARS wallThickness float 1
        //     LOOKUP_TABLE default
        //     <one value per point>
        //
        // Binary legacy VTK, XML .vtu files, vector arrays, and FIELD arrays
        // are deliberately outside this small reader's scope.
        static scalarField readVtkPointScalarField
        (
            const fileName& vtkFile,
            const word& fieldName,
            const label expectedNumberOfPoints
        );


        // Check that there is exactly one value per face and that every value
        // is finite and strictly positive.  A zero distance would create a
        // zero-volume extruded cell; a negative distance would reverse the
        // requested direction and is deliberately rejected here.
        static void validate
        (
            const scalarField& faceDistance,
            const label nSurfaceFaces
        );


        // Convert face values to point values using
        //
        //                sum(area_f * distance_f)
        // distance_p =  ----------------------------
        //                       sum(area_f)
        //
        // where the sums contain the faces connected to point p.
        static scalarField areaWeightedFaceToPoint
        (
            const MeshedSurface<face>& surface,
            const scalarField& faceDistance
        );


    // Copying is disabled to make ownership of the two fields unambiguous.
    variableDistance(const variableDistance&) = delete;
    void operator=(const variableDistance&) = delete;


public:

    // Mapping

        // Map a point field from a source surface onto arbitrary target
        // points, such as wallPatch.localPoints().
        //
        // For each target point, the function finds the closest position on
        // the source surface.  The source polygon is split into triangles and
        // the value at that closest position is calculated with barycentric
        // interpolation of the three source-vertex values.
        //
        // sourcePointValues must contain one value per source.points() entry.
        // maxDistance is a safety check: mapping stops if a target point is
        // farther than this distance from the source surface.  GREAT disables
        // this check, but a geometry-scaled tolerance is safer in production.
        //
        // An indexedOctree accelerates the nearest-triangle lookup, making the
        // method suitable for large artery surfaces such as thickMap.vtk.
        static scalarField interpolateToPoints
        (
            const MeshedSurface<face>& source,
            const scalarField& sourcePointValues,
            const pointField& targetPoints,
            const scalar maxDistance = GREAT
        );


        // Convenience wrapper for the complete operation requested here:
        //
        //  1. Read fieldName from legacy ASCII VTK POINT_DATA.
        //  2. Confirm that it has one value per source-surface point.
        //  3. Interpolate it onto targetPoints with closest-triangle and
        //     barycentric interpolation.
        //
        // The same VTK file may be used to construct source beforehand with
        // MeshedSurface<face> source(vtkFile).  Keeping source as an argument
        // avoids reading the geometry twice inside this helper.
        static scalarField interpolateVtkPointField
        (
            const MeshedSurface<face>& source,
            const fileName& vtkFile,
            const word& fieldName,
            const pointField& targetPoints,
            const scalar maxDistance = GREAT
        );

    // Constructors

        // Read face distances from distanceFile, validate them against the
        // supplied surface, and calculate point distances.
        variableDistance
        (
            const MeshedSurface<face>& surface,
            const fileName& distanceFile
        );


        // Use an already available face field.  This will be useful when the
        // VTU CellData reader is added: it can pass its values directly
        // without writing and rereading an intermediate file.
        variableDistance
        (
            const MeshedSurface<face>& surface,
            const scalarField& faceDistance
        );


    // Access

        // Return the original face-based distances without copying them.
        const scalarField& faceDistance() const;


        // Return the interpolated point-based distances without copying them.
        // This is the field that variableExtrudeMesh will use to move points.
        const scalarField& pointDistance() const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
