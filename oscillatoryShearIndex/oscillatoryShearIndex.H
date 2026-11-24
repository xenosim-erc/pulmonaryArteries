/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

Class
    Foam::functionObjects::oscillatoryShearIndex

Group
    grpForcesFunctionObjects

Description
    Computes the oscillatory shear index selected wall patches.

SourceFiles
    oscillatoryShearIndex.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_oscillatoryShearIndex_H
#define Foam_functionObjects_oscillatoryShearIndex_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class oscillatoryShearIndex Declaration
\*---------------------------------------------------------------------------*/

class oscillatoryShearIndex
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Data

        //- Write the shear stress field ?
        bool writeFields_;

        //- Accumulated wall shear stress
        //  Time integral of the wall shear stress: can be positive of negative
        volVectorField accumulatedWallShearStress_;

        //- Accumulated magnitude wall shear stress
        //  Time integral of the magnitude of the wall shear stress: can only be
        //  positive
        volScalarField accumulatedMagWallShearStress_;


protected:

    // Protected Data

        //- Wall patches to process (optionally filtered by name)
        labelList patchIDs_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate the oscillatory shear index
        void calcOscillatoryShearIndex
        (
            const volSymmTensorField& Reff,
            volScalarField& oscillatoryShearIndex,
            volVectorField& accumulatedWallShearStress,
            volScalarField& accumulatedMagWallShearStress
        );


public:

    //- Runtime type information
    TypeName("oscillatoryShearIndex");


    // Constructors

        //- Construct from Time and dictionary
        oscillatoryShearIndex
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- No copy construct
        oscillatoryShearIndex(const oscillatoryShearIndex&) = delete;

        //- No copy assignment
        void operator=(const oscillatoryShearIndex&) = delete;


    //- Destructor
    virtual ~oscillatoryShearIndex() = default;


    // Member Functions

        //- Read the oscillatoryShearIndex data
        virtual bool read(const dictionary&);

        //- Calculate the oscillatory shear index
        virtual bool execute();

        //- Report min/max and log to file,
        //- write the wall oscillatory shear index volume field.
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
