/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::PIDangularDisplacementPointPatchVectorField

Description
    Foam::PIDangularDisplacementPointPatchVectorField

SourceFiles
    PIDangularDisplacementPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef PIDangularDisplacementPointPatchVectorField_H
#define PIDangularDisplacementPointPatchVectorField_H

#include "fixedValuePointPatchField.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class PIDangularDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class PIDangularDisplacementPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private data
        //Read from pointDisplacement
        vector axis_;
        vector origin_;
        scalar angle0_;
        scalar amplitude_;
        scalar omega_;
        pointField p0_;

        // system/PIDcontrolDict
        IOdictionary PIDcontrolDict_;
        scalar P_;
        scalar I_;
        scalar D_;

        // needed vars
        scalar error_=0;
        scalar errorIntegral_=0;
        scalar oldError_=0;
        scalar oldErrorIntegral_=0;
        scalar oldOmega_=0.0;
        scalar errorDifferential_=0.0;
        label  timeIndex_=0;
        scalar auxangle_=0;
        scalar angle=0;

        int controlTarget_;
        wordList patches_;
        dictionary forcesDict_;
        scalar setPoint_=-0.35;
        label direction_;
        scalar delay_;
        List<scalar> mycoefs;
        Vector<scalar> myForce;

        //functionObjects::forceCoeffs fc;

        //bool built_=false;


public:

    //- Runtime type information
    TypeName("PIDangularDisplacement");

    // Constructors

        //- Construct from patch and internal field
        PIDangularDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        PIDangularDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        PIDangularDisplacementPointPatchVectorField
        (
            const PIDangularDisplacementPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new PIDangularDisplacementPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        PIDangularDisplacementPointPatchVectorField
        (
            const PIDangularDisplacementPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new PIDangularDisplacementPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const pointPatchFieldMapper&
            );

            //- Reverse map the given pointPatchField onto this pointPatchField
            virtual void rmap
            (
                const pointPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            IOdictionary readControl();



        //- Write
        void report(Ostream&) const;
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
