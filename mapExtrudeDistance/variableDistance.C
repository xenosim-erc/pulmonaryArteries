/*---------------------------------------------------------------------------*\
  variableDistance.C

  Definitions for the variableDistance class declared in variableDistance.H.
\*---------------------------------------------------------------------------*/

#include "variableDistance.H"

// IFstream is OpenFOAM's input-file stream.  It understands OpenFOAM list
// syntax and reports errors using the usual OpenFOAM stream machinery.
#include "IFstream.H"

// Mathematical helpers such as Foam::isfinite.
#include "scalar.H"

// OpenFOAM fatal-error reporting and exit(FatalError).
#include "error.H"

// triPointRef represents a triangle using references to three existing
// points.  Its nearestPoint() method handles projections that fall outside
// the triangle by returning the nearest point on an edge or vertex.
#include "triangle.H"

// triSurface stores triangles without changing their point labels.
// triSurfaceSearch builds an indexedOctree and answers nearest-triangle
// queries in approximately logarithmic instead of linear time.
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "pointIndexHit.H"

// std::isfinite checks that an input value is neither infinity nor NaN.
#include <cmath>

// Standard C++ streams are used for this deliberately small VTK parser.
// OpenFOAM's dictionary lexer is not appropriate for the non-OpenFOAM VTK
// header beginning with "# vtk DataFile Version ...".
#include <fstream>
#include <sstream>
#include <string>

namespace Foam
{

// * * * * * * * * * * * * Private Static Functions  * * * * * * * * * * //

scalarField variableDistance::readFaceDistance
(
    const fileName& distanceFile
)
{
    // Expand OpenFOAM/environment syntax such as "$FOAM_CASE/constant/..."
    // before attempting to open the file.
    fileName expandedFile(distanceFile);
    expandedFile.expand();

    IFstream input(expandedFile);

    if (!input.good())
    {
        FatalErrorInFunction
            << "Cannot open face-distance file " << expandedFile << nl
            << "Expected a raw OpenFOAM scalar list, for example:" << nl
            << "4" << nl
            << "(" << nl
            << "    0.001" << nl
            << "    0.0015" << nl
            << "    0.002" << nl
            << "    0.0025" << nl
            << ")" << nl
            << exit(FatalError);
    }

    // Field<T> has a stream constructor.  For scalarField it reads the list
    // length followed by the parenthesised scalar entries shown above.
    scalarField values(input);

    if (!input.good())
    {
        FatalErrorInFunction
            << "Failed while reading scalar values from " << expandedFile
            << exit(FatalError);
    }

    return values;
}


scalarField variableDistance::readVtkPointScalarField
(
    const fileName& vtkFile,
    const word& fieldName,
    const label expectedNumberOfPoints
)
{
    fileName expandedFile(vtkFile);
    expandedFile.expand();

    std::ifstream input(expandedFile.c_str());

    if (!input.good())
    {
        FatalErrorInFunction
            << "Cannot open VTK file " << expandedFile
            << exit(FatalError);
    }

    // A legacy VTK file begins with four text lines:
    //
    //   # vtk DataFile Version 3.0
    //   an arbitrary title
    //   ASCII
    //   DATASET POLYDATA
    //
    // Read these as complete lines since the title can contain spaces.
    std::string versionLine;
    std::string titleLine;
    std::string formatLine;
    std::string datasetLine;

    if
    (
        !std::getline(input, versionLine)
     || !std::getline(input, titleLine)
     || !std::getline(input, formatLine)
     || !std::getline(input, datasetLine)
    )
    {
        FatalErrorInFunction
            << expandedFile << " does not contain a complete legacy VTK "
            << "header."
            << exit(FatalError);
    }

    if (versionLine.find("vtk DataFile Version") == std::string::npos)
    {
        FatalErrorInFunction
            << expandedFile << " is not recognized as a legacy VTK file." << nl
            << "Its first line is: " << versionLine
            << exit(FatalError);
    }

    if (formatLine.find("ASCII") == std::string::npos)
    {
        FatalErrorInFunction
            << "Only legacy ASCII VTK is supported.  The format line in "
            << expandedFile << " is: " << formatLine << nl
            << "Convert binary VTK to legacy ASCII before mapping."
            << exit(FatalError);
    }

    // Search the geometry section for POINT_DATA.  The use of tokens here is
    // intentional: it lets us pass over POINTS, POLYGONS and their numeric
    // contents without needing to duplicate OpenFOAM's geometry reader.
    std::string keyword;
    label numberOfPointValues = -1;

    while (input >> keyword)
    {
        if (keyword == "POINT_DATA")
        {
            input >> numberOfPointValues;
            break;
        }
    }

    if (numberOfPointValues < 0)
    {
        FatalErrorInFunction
            << "No POINT_DATA section was found in " << expandedFile
            << exit(FatalError);
    }

    if (numberOfPointValues != expectedNumberOfPoints)
    {
        FatalErrorInFunction
            << "VTK POINT_DATA declares " << numberOfPointValues
            << " values, but the imported source surface contains "
            << expectedNumberOfPoints << " points." << nl
            << "The surface reader may have merged or reordered points, or "
            << "the field may belong to a different geometry."
            << exit(FatalError);
    }

    // A POINT_DATA section can contain multiple arrays.  Continue until the
    // requested SCALARS declaration is found.  Unknown tokens (including data
    // belonging to arrays before the requested one) are simply passed over.
    while (input >> keyword)
    {
        if (keyword != "SCALARS")
        {
            // CELL_DATA marks the end of the point-associated arrays.
            if (keyword == "CELL_DATA")
            {
                break;
            }

            continue;
        }

        std::string arrayName;
        std::string numericType;
        std::string declarationRemainder;

        input >> arrayName >> numericType;

        // The optional component count is the remainder of the SCALARS line.
        // VTK defines a default of one component when it is omitted.
        std::getline(input, declarationRemainder);
        std::istringstream declaration(declarationRemainder);

        label numberOfComponents = 1;
        declaration >> numberOfComponents;

        std::string lookupKeyword;
        std::string lookupName;
        input >> lookupKeyword >> lookupName;

        if (lookupKeyword != "LOOKUP_TABLE")
        {
            FatalErrorInFunction
                << "Expected LOOKUP_TABLE after SCALARS " << arrayName
                << " in " << expandedFile << ", but found "
                << lookupKeyword
                << exit(FatalError);
        }

        const label numberOfEntries =
            numberOfPointValues*numberOfComponents;

        if (arrayName == fieldName)
        {
            if (numberOfComponents != 1)
            {
                FatalErrorInFunction
                    << "VTK array " << fieldName << " has "
                    << numberOfComponents << " components.  A scalar "
                    << "extrusion-thickness field must have one component."
                    << exit(FatalError);
            }

            scalarField values(numberOfPointValues, Zero);

            forAll(values, pointi)
            {
                if (!(input >> values[pointi]))
                {
                    FatalErrorInFunction
                        << "Failed while reading value " << pointi
                        << " of VTK POINT_DATA array " << fieldName
                        << " from " << expandedFile
                        << exit(FatalError);
                }

                if (!std::isfinite(values[pointi]))
                {
                    FatalErrorInFunction
                        << "VTK POINT_DATA array " << fieldName
                        << " contains a non-finite value at point " << pointi
                        << ": " << values[pointi]
                        << exit(FatalError);
                }
            }

            Info<< "Read VTK POINT_DATA scalar array " << fieldName
                << " with " << values.size() << " values from "
                << expandedFile << endl;

            return values;
        }

        // This is a different scalar array.  Consume all of its values before
        // resuming the search for another SCALARS declaration.
        scalar ignoredValue = Zero;

        for (label entryi = 0; entryi < numberOfEntries; ++entryi)
        {
            if (!(input >> ignoredValue))
            {
                FatalErrorInFunction
                    << "Failed while skipping VTK scalar array " << arrayName
                    << " in " << expandedFile
                    << exit(FatalError);
            }
        }
    }

    FatalErrorInFunction
        << "Could not find scalar POINT_DATA array " << fieldName
        << " in " << expandedFile << nl
        << "Expected a declaration such as:" << nl
        << "SCALARS " << fieldName << " float 1" << nl
        << "LOOKUP_TABLE default"
        << exit(FatalError);

    // Unreachable, but required to satisfy the C++ return type.
    return scalarField();
}


void variableDistance::validate
(
    const scalarField& faceDistance,
    const label nSurfaceFaces
)
{
    // A positional face field is meaningful only if its size exactly matches
    // the number of faces in the surface.
    if (faceDistance.size() != nSurfaceFaces)
    {
        FatalErrorInFunction
            << "The distance field contains " << faceDistance.size()
            << " values, but the surface contains " << nSurfaceFaces
            << " faces." << nl
            << "There must be exactly one distance value per surface face,"
            << " in identical face order."
            << exit(FatalError);
    }

    forAll(faceDistance, facei)
    {
        const scalar distance = faceDistance[facei];

        if (!std::isfinite(distance))
        {
            FatalErrorInFunction
                << "Distance for surface face " << facei
                << " is not finite: " << distance
                << exit(FatalError);
        }

        if (distance <= 0)
        {
            FatalErrorInFunction
                << "Distance for surface face " << facei
                << " is not positive: " << distance << nl
                << "Positive distances are required to avoid zero-volume or "
                << "reversed cells. Reverse the surface normals (or use "
                << "flipNormals) to change extrusion direction."
                << exit(FatalError);
        }
    }
}


scalarField variableDistance::areaWeightedFaceToPoint
(
    const MeshedSurface<face>& surface,
    const scalarField& faceDistance
)
{
    // Both accumulators have one entry per surface point and start at zero.
    scalarField weightedDistance(surface.points().size(), Zero);
    scalarField totalArea(surface.points().size(), Zero);

    const List<face>& faces = surface.surfFaces();

    // MeshedSurface caches the magnitude of every face area.  Area weighting
    // makes a tiny neighbouring face influence a shared point less than a
    // large face.  For equal-sized faces this reduces to an ordinary average.
    const scalarField& faceArea = surface.magSf();

    forAll(faces, facei)
    {
        const face& currentFace = faces[facei];
        const scalar area = faceArea[facei];

        if (area <= VSMALL)
        {
            FatalErrorInFunction
                << "Surface face " << facei
                << " has zero or near-zero area (" << area << ")." << nl
                << "Degenerate faces must be repaired before extrusion."
                << exit(FatalError);
        }

        // A face stores the global surface-point labels at its vertices.
        // Add this face's contribution to each one of those points.
        forAll(currentFace, facePointi)
        {
            const label pointi = currentFace[facePointi];

            weightedDistance[pointi] += area*faceDistance[facei];
            totalArea[pointi] += area;
        }
    }

    forAll(weightedDistance, pointi)
    {
        if (totalArea[pointi] <= VSMALL)
        {
            FatalErrorInFunction
                << "Surface point " << pointi
                << " is not connected to any non-degenerate face."
                << exit(FatalError);
        }

        weightedDistance[pointi] /= totalArea[pointi];
    }

    // The accumulator now contains the completed point-distance field.
    return weightedDistance;
}


// * * * * * * * * * * * Public Mapping Function  * * * * * * * * * * * //

scalarField variableDistance::interpolateToPoints
(
    const MeshedSurface<face>& source,
    const scalarField& sourcePointValues,
    const pointField& targetPoints,
    const scalar maxDistance
)
{
    const pointField& sourcePoints = source.points();
    const List<face>& sourceFaces = source.surfFaces();

    if (sourcePointValues.size() != sourcePoints.size())
    {
        FatalErrorInFunction
            << "The source point field contains "
            << sourcePointValues.size() << " values, but the source surface "
            << "contains " << sourcePoints.size() << " points." << nl
            << "A VTK POINT_DATA field must contain exactly one value per "
            << "VTK POINTS entry."
            << exit(FatalError);
    }

    if (sourceFaces.empty())
    {
        FatalErrorInFunction
            << "Cannot interpolate from a surface with no faces."
            << exit(FatalError);
    }

    if (maxDistance <= 0)
    {
        FatalErrorInFunction
            << "maxDistance must be positive; received " << maxDistance
            << exit(FatalError);
    }

    // Convert the possibly polygonal MeshedSurface to triangles.  Each fan
    // triangle retains the original source point labels, which means a hit
    // triangle can index sourcePointValues directly without a separate map.
    label numberOfTriangles = 0;

    forAll(sourceFaces, facei)
    {
        if (sourceFaces[facei].size() < 3)
        {
            FatalErrorInFunction
                << "Source face " << facei << " has only "
                << sourceFaces[facei].size() << " vertices."
                << exit(FatalError);
        }

        numberOfTriangles += sourceFaces[facei].size() - 2;
    }

    triFaceList triangles(numberOfTriangles);
    label trianglei = 0;

    forAll(sourceFaces, facei)
    {
        const face& sourceFace = sourceFaces[facei];

        for (label fp = 1; fp + 1 < sourceFace.size(); ++fp)
        {
            triangles[trianglei++] = triFace
            (
                sourceFace[0],
                sourceFace[fp],
                sourceFace[fp + 1]
            );
        }
    }

    const triSurface searchSurface(triangles, sourcePoints);
    const triSurfaceSearch search(searchSurface);

    // The octree is constructed on first use.  Building it once here avoids
    // the previous O(numberOfTargets * numberOfTriangles) exhaustive search.
    const indexedOctree<treeDataTriSurface>& tree = search.tree();

    scalarField targetValues(targetPoints.size(), Zero);
    scalar largestMappingDistance = 0;

    // findNearest expects a squared search radius.  GREAT means effectively
    // unbounded; otherwise the caller's physical tolerance is squared.
    const scalar searchDistanceSqr =
        (maxDistance >= sqrt(GREAT) ? GREAT : sqr(maxDistance));

    forAll(targetPoints, targetPointi)
    {
        const point& target = targetPoints[targetPointi];
        const pointIndexHit hit = tree.findNearest(target, searchDistanceSqr);

        if (!hit.hit())
        {
            FatalErrorInFunction
                << "No source triangle was found within maxDistance "
                << maxDistance << " of target point " << targetPointi
                << " at " << target << nl
                << "Check that thickMap.vtk and the OpenFOAM wall patch use "
                << "the same coordinates and units."
                << exit(FatalError);
        }

        const triFace& triangle = searchSurface[hit.index()];
        const point& a = sourcePoints[triangle[0]];
        const point& b = sourcePoints[triangle[1]];
        const point& c = sourcePoints[triangle[2]];
        const point& closest = hit.point();

        const vector edge0 = b - a;
        const vector edge1 = c - a;
        const vector relative = closest - a;
        const scalar d00 = edge0 & edge0;
        const scalar d01 = edge0 & edge1;
        const scalar d11 = edge1 & edge1;
        const scalar denominator = d00*d11 - d01*d01;

        if (denominator <= VSMALL)
        {
            FatalErrorInFunction
                << "Octree returned zero-area source triangle "
                << hit.index() << "."
                << exit(FatalError);
        }

        scalar weight1 =
            (d11*(relative & edge0) - d01*(relative & edge1))/denominator;
        scalar weight2 =
            (d00*(relative & edge1) - d01*(relative & edge0))/denominator;
        scalar weight0 = 1 - weight1 - weight2;

        // Protect interpolation on triangle edges from floating-point drift.
        weight0 = max(scalar(0), min(scalar(1), weight0));
        weight1 = max(scalar(0), min(scalar(1), weight1));
        weight2 = max(scalar(0), min(scalar(1), weight2));

        const scalar weightSum = weight0 + weight1 + weight2;
        weight0 /= weightSum;
        weight1 /= weightSum;
        weight2 /= weightSum;

        targetValues[targetPointi] =
            weight0*sourcePointValues[triangle[0]]
          + weight1*sourcePointValues[triangle[1]]
          + weight2*sourcePointValues[triangle[2]];

        largestMappingDistance = max
        (
            largestMappingDistance,
            mag(target - closest)
        );
    }

    Info<< "Mapped " << targetPoints.size()
        << " target points from the source surface." << nl
        << "Largest source-surface projection distance: "
        << largestMappingDistance << endl;

    return targetValues;
}


scalarField variableDistance::interpolateVtkPointField
(
    const MeshedSurface<face>& source,
    const fileName& vtkFile,
    const word& fieldName,
    const pointField& targetPoints,
    const scalar maxDistance
)
{
    const scalarField sourcePointValues
    (
        readVtkPointScalarField
        (
            vtkFile,
            fieldName,
            source.points().size()
        )
    );

    return interpolateToPoints
    (
        source,
        sourcePointValues,
        targetPoints,
        maxDistance
    );
}


// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

variableDistance::variableDistance
(
    const MeshedSurface<face>& surface,
    const fileName& distanceFile
)
:
    // Members are initialized in the order in which they are declared in the
    // header.  First read the face field, then use it to build the point field.
    faceDistance_(readFaceDistance(distanceFile)),
    pointDistance_()
{
    validate(faceDistance_, surface.size());
    pointDistance_ = areaWeightedFaceToPoint(surface, faceDistance_);
}


variableDistance::variableDistance
(
    const MeshedSurface<face>& surface,
    const scalarField& faceDistance
)
:
    // scalarField owns its data, so this makes an independent copy of the
    // caller's field.  The class therefore remains valid if the caller's
    // temporary VTU data is later destroyed.
    faceDistance_(faceDistance),
    pointDistance_()
{
    validate(faceDistance_, surface.size());
    pointDistance_ = areaWeightedFaceToPoint(surface, faceDistance_);
}


// * * * * * * * * * * * * * Access Functions * * * * * * * * * * * * * //

const scalarField& variableDistance::faceDistance() const
{
    return faceDistance_;
}


const scalarField& variableDistance::pointDistance() const
{
    return pointDistance_;
}

} // End namespace Foam

// ************************************************************************* //
