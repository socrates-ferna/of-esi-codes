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
#include "lforces.H"
#include "lforceCoeffs.H"


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
    angle0_(Zero),
    amplitude_(Zero),
    omega_(Zero),
    p0_(p.localPoints()),
    PIDcontrolDict_(this->readControl()),
    P_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("P",0.5)),
    I_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("I",0.5)),
    D_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("D",1.0)),
    controlTarget_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault<int>("controlTarget",2)),
    forcesDict_(PIDcontrolDict_.subDict("PIDcontroller").subDict("controlledVarData")),
    setPoint_(forcesDict_.getOrDefault<scalar>("setPoint",0.35)),
    direction_(forcesDict_.getOrDefault<scalar>("direction",3))
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
    D_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("D",1.0)),
    controlTarget_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault<int>("controlTarget",2)),
    forcesDict_(PIDcontrolDict_.subDict("PIDcontroller").subDict("controlledVarData")),
    setPoint_(forcesDict_.getOrDefault<scalar>("setPoint",0.35)),
    direction_(forcesDict_.getOrDefault<label>("direction",3))
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
    PIDcontrolDict_(ptf.PIDcontrolDict_),
    P_(ptf.P_),
    I_(ptf.I_),
    D_(ptf.D_),
    controlTarget_(ptf.controlTarget_),
    forcesDict_(ptf.forcesDict_),
    setPoint_(ptf.setPoint_),
    direction_(ptf.direction)
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
    PIDcontrolDict_(ptf.PIDcontrolDict_),
    P_(ptf.P_),
    I_(ptf.I_),
    D_(ptf.D_),
    controlTarget_(ptf.controlTarget_),
    forcesDict_(ptf.forcesDict_),
    setPoint_(ptf.setPoint_)
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

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    if (timeIndex_!= t.timeIndex())
    {
        timeIndex_ = t.timeIndex();
        oldOmega_= omega_;
        oldError_ = error_;
        oldErrorIntegral_ = errorIntegral_;
    }

    functionObjects::lforces f("forces",t,forcesDict_,true);
    Vector<float> myForce;  // vector before. Vector is OF 
    List<scalar> mycoefs(6);
    /*if (controlTarget_ == 0 || controlTarget_ == 1)
    {
        Info<<"   ";
        functionObjects::forceCoeffs fc("forceCoeffs",t,forcesDict_,true);
    }*/

    /*
    functionObjects::forceCoeffs fc = 
    functionObjects::forceCoeffs ("forceCoeffs",t,forcesDict_,true);
    */

   functionObjects::lforceCoeffs fc("forceCoeffs",t,forcesDict_,true);

    switch (controlTarget_) {
        case 0:
            // not implemented yet
            
            //functionObjects::forces f("forces", mesh, forcesDict_);
            f.calcForcesMoment();

            myForce = f.forceEff();
            error_ = setPoint_ - myForce.x(); // this should choose x,y,z according to direction

            errorIntegral_ = oldErrorIntegral_ + I_*0.5*(error_ + oldError_)*t.deltaTValue();
            errorDifferential_ = oldError_ - error_;
            //const scalar errorDifferential = oldError_ - error_;
            omega_ = oldOmega_ + P_*error_ + errorIntegral_ + D_*errorDifferential_/t.deltaTValue();


            Info<< "totalForce is: "<< myForce<<endl;
            Info<< "dir1 force is: "<< myForce.x() <<endl;
            /*
            Info<< "comp0 force is: "<< myForce.component(0) <<endl;
            Info<< "comp1 force is: "<< myForce.component(1) <<endl;
            Info<< "comp2 force is: "<< myForce.component(2) <<endl;
            */
            Info<< "target force is: "<< setPoint_ <<endl;
            //Info<< "e1:"<< f.coordSys_.e1()<<endl; NO VALE PORQUE SON PRIVADAS
            //Info<< "e2:"<< f.coordSys_.e2()<<endl;
            //Info<< "e3:"<< f.coordSys_.e3()<<endl;

        break;

        case 1:  // PUEDO ACCEDER AL OBJETO YA CREADO POR PIMPLE?
            // not implemented yet

            // IMPLEMENT YOUR OWN FORCECOEFFS OR SIMPLY CALL EXECUTE?
            Info<<"forceCoeffs"<<endl;
            //fc.calcForcesMoment();
            fc.calcrot();
            mycoefs = fc.coefList;
            Info<<"Cd"<< mycoefs[0]<<endl;
            Info<<"Cs"<< mycoefs[1]<<endl;
            Info<<"Cl"<< mycoefs[2]<<endl;
            Info<<"CmRoll"<< mycoefs[3]<<endl;
            Info<<"CmPitch"<< mycoefs[4]<<endl;
            Info<<"CmYaw"<< mycoefs[5]<<endl;

            //fc.execute();
            /*
            Info<< "COEFStotalForce is: "<< myForce<<endl;
            Info<< "COEFSdir1 force is: "<< myForce.x() <<endl;
            Info<< "target force is: "<< setPoint_ <<endl;
            */
            /*
            f.calcForcesMoment();
            
            List<Field<scalar>> liftCoeffs(3);
            scalar ClTot = 0;

            const scalar pDyn = f.rhoRef_*sqr(f.maUInf_); // this won't work, private data members

            const scalar forceScaling = 1.0/(Aref_*pDyn + SMALL);

            forAll(liftCoeffs, i)
            {
                const Field<vector> localForce(f.coordSys_.localVector(myForce[i]));
            }
            */
            break;

        case 2:

            Info<< "angle antes: "<< angle<< endl;
            //Info<< dMD.lookup("dummyValue")<< endl;
            //write(Info);

            // PID CONTROL
            error_ = setPoint_ - angle;

            errorIntegral_ = oldErrorIntegral_ + I_*0.5*(error_ + oldError_)*t.deltaTValue();
            //const scalar errorDifferential = oldError_ - error_;
            errorDifferential_ = oldError_ - error_;
            omega_ = oldOmega_ + P_*error_ + errorIntegral_ + D_*errorDifferential_/t.deltaTValue();
            break;

        default:

            Info<< "Unknown controlled variable"<<endl
            <<"Allowed variables are:"<< endl<< "(forces forceCoeffs angle)"
            <<endl<<exit(FatalError);

        break;

        // EXTERNAL CONTROL WILL BE DONE IN ANOTHER pointPatchFields option
    }


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
