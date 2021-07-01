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

\*---------------------------------------------------------------------------*/

#include "angularDisplacementPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

angularDisplacementPointPatchVectorField::
angularDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    axis_(Zero),
    origin_(Zero),
    angle0_(0.0),
    amplitude_(0.0),
    omega_(0.0),
    p0_(p.localPoints())
{}


angularDisplacementPointPatchVectorField::
angularDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    axis_(dict.lookup("axis")),
    origin_(dict.lookup("origin")),
    angle0_(dict.get<scalar>("angle0")),
    amplitude_(dict.get<scalar>("amplitude")),
    omega_(dict.get<scalar>("omega"))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("p0"))
    {
        p0_ = vectorField("p0", dict , p.size());
    }
    else
    {
        p0_ = p.localPoints();
    }
}


angularDisplacementPointPatchVectorField::
angularDisplacementPointPatchVectorField
(
    const angularDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    axis_(ptf.axis_),
    origin_(ptf.origin_),
    angle0_(ptf.angle0_),
    amplitude_(ptf.amplitude_),
    omega_(ptf.omega_),
    p0_(ptf.p0_, mapper)
{}


angularDisplacementPointPatchVectorField::
angularDisplacementPointPatchVectorField
(
    const angularDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    axis_(ptf.axis_),
    origin_(ptf.origin_),
    angle0_(ptf.angle0_),
    amplitude_(ptf.amplitude_),
    omega_(ptf.omega_),
    p0_(ptf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void angularDisplacementPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);

    p0_.autoMap(m);
}


void angularDisplacementPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const angularDisplacementPointPatchVectorField& aODptf =
        refCast<const angularDisplacementPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(aODptf, addr);

    p0_.rmap(aODptf.p0_, addr);
}


void angularDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    ++testlabel;
    Info<< "el contador vale:"<< testlabel<<endl;
    //this->write(Info);
    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    // LEE EL DYNAMICMESHDICT IFMODIFIED IOdictionary dMD = mesh.lookupObject<IOdictionary>("dynamicMeshDict");
    Info<< "angle antes: "<< angle<< endl;
    //Info<< dMD.lookup("dummyValue")<< endl;
    angle = angle0_ + omega_*t.value();  //amplitude_*sin(omega_*t.value());
    vector axisHat = axis_/mag(axis_);
    vectorField p0Rel(p0_ - origin_);
    Info<< "angle después"<< angle<< endl;
    vectorField::operator=
    (
        p0Rel*(cos(angle) - 1)
      + (axisHat ^ p0Rel*sin(angle))
      + (axisHat & p0Rel)*(1 - cos(angle))*axisHat
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
    //Info<< "db access"<< db().names() << endl;

}


void angularDisplacementPointPatchVectorField::write
(
    Ostream& os
) const
{
    pointPatchField<vector>::write(os);
    os.writeEntry("axis", axis_);
    os.writeEntry("origin", origin_);
    os.writeEntry("angle0", angle0_);
    os.writeEntry("angle", angle);
    os.writeEntry("amplitude", amplitude_);
    os.writeEntry("omega", omega_);
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    angularDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
