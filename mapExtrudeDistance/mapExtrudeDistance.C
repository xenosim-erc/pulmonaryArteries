/*---------------------------------------------------------------------------*\
  mapExtrudeDistance

  Read a scalar POINT_DATA array from a legacy ASCII VTK surface, interpolate
  it onto an OpenFOAM boundary patch, and write it as a pointScalarField.

  This utility changes no mesh geometry and performs no extrusion.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "MeshedSurfaces.H"
#include "polyPatch.H"
#include "variableDistance.H"

using namespace Foam;

int main(int argc, char* argv[])
{
    argList::addNote
    (
        "Map legacy-ASCII VTK POINT_DATA onto an OpenFOAM boundary patch "
        "and write a pointScalarField. No mesh points are moved."
    );

    argList::addOption
    (
        "vtk",
        "file",
        "VTK surface file (default: thickMap.vtk in the case directory)"
    );
    argList::addOption
    (
        "patch",
        "name",
        "Target OpenFOAM boundary patch (default: wall)"
    );
    argList::addOption
    (
        "vtkField",
        "name",
        "VTK POINT_DATA scalar array (default: WallThickness)"
    );
    argList::addOption
    (
        "outputField",
        "name",
        "OpenFOAM pointScalarField name (default: WallThickness)"
    );
    argList::addOption
    (
        "maxDistance",
        "scalar",
        "Maximum VTK-to-patch projection distance (default: 0.2)"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    fileName vtkFile
    (
        args.getOrDefault<fileName>("vtk", runTime.path()/"thickMap.vtk")
    );

    // Interpret relative VTK names relative to the selected -case directory,
    // not relative to the directory from which this executable was launched.
    if (!vtkFile.isAbsolute())
    {
        vtkFile = runTime.path()/vtkFile;
    }
    const word patchName(args.getOrDefault<word>("patch", "wall"));
    const word vtkFieldName
    (
        args.getOrDefault<word>("vtkField", "WallThickness")
    );
    const word outputFieldName
    (
        args.getOrDefault<word>("outputField", vtkFieldName)
    );
    const scalar maxDistance
    (
        args.getOrDefault<scalar>("maxDistance", 0.2)
    );

    const label patchi = mesh.boundaryMesh().findPatchID(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Cannot find patch " << patchName << "." << nl
            << "Available patches are: " << mesh.boundaryMesh().names()
            << exit(FatalError);
    }

    const polyPatch& targetPatch = mesh.boundaryMesh()[patchi];

    Info<< "Reading VTK surface " << vtkFile << nl
        << "Mapping POINT_DATA array " << vtkFieldName << nl
        << "Target patch " << patchName << " has "
        << targetPatch.nPoints() << " points and "
        << targetPatch.size() << " faces" << nl
        << "Maximum mapping distance: " << maxDistance << nl << endl;

    // MeshedSurface reads the VTK points and polygon connectivity.  The
    // variableDistance helper separately parses the named POINT_DATA array.
    const MeshedSurface<face> vtkSurface(vtkFile);

    const scalarField patchValues
    (
        variableDistance::interpolateVtkPointField
        (
            vtkSurface,
            vtkFile,
            vtkFieldName,
            targetPatch.localPoints(),
            maxDistance
        )
    );

    const pointMesh& pMesh = pointMesh::New(mesh);

    // A pointScalarField has one primitive value for every global mesh point.
    // Values not belonging to the selected wall patch remain zero.  The field
    // carries length dimensions because WallThickness is a distance.
    pointScalarField mappedField
    (
        IOobject
        (
            outputFieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh,
        dimensionedScalar("zero", dimLength, Zero)
    );

    scalarField& allPointValues = mappedField.primitiveFieldRef();
    const labelList& patchMeshPoints = targetPatch.meshPoints();

    forAll(patchMeshPoints, patchPointi)
    {
        allPointValues[patchMeshPoints[patchPointi]] = patchValues[patchPointi];
    }

    // Populate point-patch representations from the global point values.
    mappedField.correctBoundaryConditions();

    if (!mappedField.write())
    {
        FatalErrorInFunction
            << "Failed to write point field " << outputFieldName
            << exit(FatalError);
    }

    Info<< nl << "Mapped field statistics on patch " << patchName << ':' << nl
        << "    minimum = " << min(patchValues) << nl
        << "    maximum = " << max(patchValues) << nl
        << "    average = " << average(patchValues) << nl
        << "Wrote " << mappedField.objectPath() << nl
        << "Mesh geometry was not modified." << nl << endl;

    return 0;
}

// ************************************************************************* //
