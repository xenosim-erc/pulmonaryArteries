/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
    Modified 2026-11-24 for the pulmonaryArteries repository: the
    accumulated (time-integrated) WSS and |WSS| fields are read back when the
    solver restarts (READ_IF_PRESENT instead of NO_READ). Taken from
    https://github.com/solids4foam/oscillatoryShearIndex (commit 57cd043).
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "oscillatoryShearIndex.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(oscillatoryShearIndex, 0);
    addToRunTimeSelectionTable(functionObject, oscillatoryShearIndex, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::oscillatoryShearIndex::writeFileHeader
(
    Ostream& os
) const
{
    writeHeader(os, "Oscillatory shear index");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    os  << endl;
}


void Foam::functionObjects::oscillatoryShearIndex::calcOscillatoryShearIndex
(
    const volSymmTensorField& Reff,
    volScalarField& oscillatoryShearIndex,
    volVectorField& accumulatedWallShearStress,
    volScalarField& accumulatedMagWallShearStress
)
{
    // Time step
    const scalar deltaT = Reff.time().deltaTValue();

    for (const label patchi : patchIDs_)
    {
        scalarField& osi = oscillatoryShearIndex.boundaryFieldRef()[patchi];
        vectorField& accumWSS =
            accumulatedWallShearStress.boundaryFieldRef()[patchi];
        const vectorField& accumWSSOld =
            accumulatedWallShearStress.oldTime().boundaryFieldRef()[patchi];
        scalarField& accumMagWSS =
            accumulatedMagWallShearStress.boundaryFieldRef()[patchi];
        const scalarField& accumMagWSSOld =
            accumulatedMagWallShearStress.oldTime().boundaryFieldRef()[patchi];
        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];
        const symmTensorField& Reffp = Reff.boundaryField()[patchi];

        // Current wall shear stress
        const vectorField wss((-Sfp/magSfp) & Reffp);

        // Wall shear stress increment
        const vectorField deltaWSS(wss*deltaT);

        // Update the accumulated wall shear stress
        accumWSS = accumWSSOld + deltaWSS;

        // Update the accumulated magnitude wall shear stress
        accumMagWSS = accumMagWSSOld + mag(deltaWSS);

        // Calculate the oscillatory shear index
        osi = 0.5*(1.0 - mag(accumWSS)/accumMagWSS);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::oscillatoryShearIndex::oscillatoryShearIndex
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    writeFields_(true),  // May change in the future
    accumulatedWallShearStress_
    (
        IOobject
        (
            scopedName(typeName) + "_accumulatedWallShearStress",
            mesh_.time().timeName(),
            mesh_,
            // Read back on a restart, so the integrals carry on from the
            // start time instead of starting again from zero
            IOobjectOption::READ_IF_PRESENT,
            IOobjectOption::NO_WRITE
        ),
        mesh_,
        dimensionedVector(sqr(dimLength)/dimTime, Zero)
    ),
    accumulatedMagWallShearStress_
    (
        IOobject
        (
            scopedName(typeName) + "_accumulatedMagWallShearStress",
            mesh_.time().timeName(),
            mesh_,
            // Read back on a restart, so the integrals carry on from the
            // start time instead of starting again from zero
            IOobjectOption::READ_IF_PRESENT,
            IOobjectOption::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(sqr(dimLength)/dimTime, 0.0)
    )
{
    read(dict);

    writeFileHeader(file());

    volScalarField* oscillatoryShearIndexPtr
    (
        new volScalarField
        (
            IOobject
            (
                scopedName(typeName),
                mesh_.time().timeName(),
                mesh_,
                IOobjectOption::NO_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::REGISTER
            ),
            mesh_,
            dimensionedScalar(sqr(dimLength)/sqr(dimTime), Zero)
        )
    );

    mesh_.objectRegistry::store(oscillatoryShearIndexPtr);

    // Store old time fields
    accumulatedWallShearStress_.storeOldTime();
    accumulatedMagWallShearStress_.storeOldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::oscillatoryShearIndex::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    writeFields_ = true;   // May change in the future
    dict.readIfPresent("writeFields", writeFields_);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    wordRes patchNames;
    labelHashSet patchSet(0);
    if (dict.readIfPresent("patches", patchNames) && !patchNames.empty())
    {
        patchSet = pbm.patchSet(patchNames);
    }

    labelHashSet allWalls(pbm.findPatchIDs<wallPolyPatch>());

    Info<< type() << ' ' << name() << ':' << nl;

    if (patchSet.empty())
    {
        patchIDs_ = allWalls.sortedToc();

        Info<< "    processing all (" << patchIDs_.size()
            << ") wall patches" << nl << endl;
    }
    else
    {
        allWalls &= patchSet;
        patchSet -= allWalls;
        patchIDs_ = allWalls.sortedToc();

        if (!patchSet.empty())
        {
            WarningInFunction
                << "Requested wall shear stress on ("
                << patchSet.size() << ") non-wall patches:" << nl;

            for (const label patchi : patchSet.sortedToc())
            {
                Info<< "        " << pbm[patchi].name() << nl;
            }
            Info<< nl;
        }

        Info<< "    processing (" << patchIDs_.size()
            << ") wall patches:" << nl;

        for (const label patchi : patchIDs_)
        {
            Info<< "        " << pbm[patchi].name() << nl;
        }
        Info<< endl;
    }

    return true;
}


bool Foam::functionObjects::oscillatoryShearIndex::execute()
{
    auto& oscillatoryShearIndex =
        mesh_.lookupObjectRef<volScalarField>(scopedName(typeName));

    // Compressible
    {
        typedef compressible::turbulenceModel turbType;

        const turbType* modelPtr =
            findObject<turbType>(turbulenceModel::propertiesName);

        if (modelPtr)
        {
            calcOscillatoryShearIndex
            (
                modelPtr->devRhoReff(),
                oscillatoryShearIndex,
                accumulatedWallShearStress_,
                accumulatedMagWallShearStress_
            );
            return true;
        }
    }

    // Incompressible
    {
        typedef incompressible::turbulenceModel turbType;

        const turbType* modelPtr =
            findObject<turbType>(turbulenceModel::propertiesName);

        if (modelPtr)
        {
            calcOscillatoryShearIndex
            (
                modelPtr->devReff(),
                oscillatoryShearIndex,
                accumulatedWallShearStress_,
                accumulatedMagWallShearStress_
            );
            return true;
        }
    }

    FatalErrorInFunction
        << "Unable to find turbulence model in the "
        << "database" << exit(FatalError);

    return false;
}


bool Foam::functionObjects::oscillatoryShearIndex::write()
{
    const auto& oscillatoryShearIndex =
        obr_.lookupObject<volScalarField>(scopedName(typeName));

    Log << type() << ' ' << name() << " write:" << nl;

    if (writeFields_)
    {
        Log << "    writing field " << oscillatoryShearIndex.name() << endl;
        oscillatoryShearIndex.write();

        Log << "    writing field " << accumulatedWallShearStress_.name() << endl;
        accumulatedWallShearStress_.write();

        Log << "    writing field " << accumulatedMagWallShearStress_.name() << endl;
        accumulatedMagWallShearStress_.write();
    }

    const fvPatchList& patches = mesh_.boundary();

    for (const label patchi : patchIDs_)
    {
        const fvPatch& pp = patches[patchi];

        const scalarField& osi = oscillatoryShearIndex.boundaryField()[patchi];

        const MinMax<scalar> limits = gMinMax(osi);

        if (UPstream::master())
        {
            writeCurrentTime(file());

            file()
                << token::TAB << pp.name()
                << token::TAB << limits.min()
                << token::TAB << limits.max()
                << endl;
        }

        Log << "    min/max(" << pp.name() << ") = "
            << limits.min() << ", " << limits.max() << endl;
    }

    return true;
}


// ************************************************************************* //
