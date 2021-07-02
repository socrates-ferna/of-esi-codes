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

#include "PIDangularDisplacementPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"
#include "forces.H"
#include "forceCoeffs.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

PIDangularDisplacementPointPatchVectorField::
PIDangularDisplacementPointPatchVectorField
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
    p0_(p.localPoints()),
    PIDcontrolDict_(this->readControl()),
    P_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("P",0.5)),
    I_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("I",0.5)),
    D_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("D",1.0))
{}


PIDangularDisplacementPointPatchVectorField::
PIDangularDisplacementPointPatchVectorField
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
    omega_(dict.get<scalar>("omega")),
    PIDcontrolDict_(this->readControl()),
    P_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("P",0.5)),
    I_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("I",0.5)),
    D_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("D",1.0))
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
    //PIDcontrolDict_ = readControl();
    
}


PIDangularDisplacementPointPatchVectorField::
PIDangularDisplacementPointPatchVectorField
(
    const PIDangularDisplacementPointPatchVectorField& ptf,
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
    p0_(ptf.p0_, mapper),
    PIDcontrolDict_(this->readControl()),
    P_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("P",0.5)),
    I_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("I",0.5)),
    D_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("D",1.0))
{}


PIDangularDisplacementPointPatchVectorField::
PIDangularDisplacementPointPatchVectorField
(
    const PIDangularDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    axis_(ptf.axis_),
    origin_(ptf.origin_),
    angle0_(ptf.angle0_),
    amplitude_(ptf.amplitude_),
    omega_(ptf.omega_),
    p0_(ptf.p0_),
    PIDcontrolDict_(this->readControl()),
    P_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("P",0.5)),
    I_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("I",0.5)),
    D_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("D",1.0))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void PIDangularDisplacementPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);

    p0_.autoMap(m);
}


void PIDangularDisplacementPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const PIDangularDisplacementPointPatchVectorField& aODptf =
        refCast<const PIDangularDisplacementPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(aODptf, addr);

    p0_.rmap(aODptf.p0_, addr);
}


void PIDangularDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    /*
    if (!built_)
    {
        IOdictionary PIDcontrolDict_
        (
            IOobject
            (
                "PIDcontrolDict",
                this->internalField().mesh().time().system(),
                this->internalField().mesh()(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::AUTO_WRITE
            )
        );
        built_=true;
        
        Info<< "READ PIDCONTROLDICT:"<< endl<< PIDcontrolDict_.tokens() << endl;

        readControl();

    }
    */

   //Info<< "P: "<< PIDcontrolDict_.subDict("PIDcontroller").lookup("P")<< endl;
   //Info<< "mesh inside PointPatch, ObjectRegistry: "<< endl<< mesh.names()<< endl;

    //++testlabel;
    //Info<< "el contador vale:"<< testlabel<<endl;
    //this->write(Info);

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    if (timeIndex_!= t.timeIndex())
    {
        timeIndex_ = t.timeIndex();
        oldOmega_= omega_;
        oldError_ = error_;
        oldErrorIntegral_ = errorIntegral_;
    }
    switch (controlTarget_) {
        case "forces":

        break

        case "forceCoeffs":

        break

        case "angle":

        break
    }
    if (control)

    Info<< "angle antes: "<< angle<< endl;
    //Info<< dMD.lookup("dummyValue")<< endl;
    //write(Info);

    // PID CONTROL
    error_ = setPoint - angle;

    errorIntegral_ = oldErrorIntegral_ + I_*0.5*(error_ + oldError_)*t.deltaTValue();
    const scalar errorDifferential = oldError_ - error_;
    omega_ = oldOmega_ + P_*error_ + errorIntegral_ + D_*errorDifferential/t.deltaTValue();
    // MOVEMENT BASED ON CALCULATED OMEGA
    angle = angle0_ + angle + omega_*t.deltaTValue();
    //auxangle_ = angle0_ + omega_*t.value();  //amplitude_*sin(omega_*t.value());
    vector axisHat = axis_/mag(axis_);
    vectorField p0Rel(p0_ - origin_);
    // report calculated angles
    report(Info);
    vectorField::operator=
    (
        p0Rel*(cos(angle) - 1)
      + (axisHat ^ p0Rel*sin(angle))
      + (axisHat & p0Rel)*(1 - cos(angle))*axisHat
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
    //Info<< "db access"<< db().names() << endl;

}


void PIDangularDisplacementPointPatchVectorField::write
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
    //p0_.writeEntry("p0", os);
    os.writeEntry("P",P_);
    os.writeEntry("I",I_);
    os.writeEntry("D",D_);
    writeEntry("value", os);
}

IOdictionary PIDangularDisplacementPointPatchVectorField::readControl()
{
    IOdictionary readdict
    (
        IOobject
        (
            "PIDcontrolDict",
            this->internalField().mesh().time().system(),
            this->internalField().mesh()(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::AUTO_WRITE
        )
    );
    
    Info<< "READCONTROLS FUNCTION PIDCONTROLDICT:"<< endl<< readdict.dictName() << endl;
    return readdict;
}

void PIDangularDisplacementPointPatchVectorField::report
(
    Ostream& os
) const
{
    os.writeEntry("oldOmega", oldOmega_);
    os.writeEntry("omega", omega_);
    os.writeEntry("error", error_);
    os.writeEntry("errorIntegral", errorIntegral_);
    //os.writeEntry("auxangle", auxangle_);
    os.writeEntry("angle", angle);
    return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    PIDangularDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
