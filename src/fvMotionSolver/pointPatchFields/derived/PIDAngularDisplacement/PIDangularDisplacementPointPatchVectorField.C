/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
        Author: Socrates Fernandez
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
#include "FIFOStack.H" // implementation pending

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
    controlDelay_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("controlDelay",1.0)),
    controlTarget_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault<int>("controlTarget",2)),
    forcesDict_(PIDcontrolDict_.subDict("PIDcontroller").subDict("controlledVarData")),
    setPoint_(forcesDict_.getOrDefault<scalar>("setPoint",0.35)),
    direction_(forcesDict_.getOrDefault<scalar>("direction",3)),
    mycoefs(6,Zero),
    myForce(Zero),
    anglemax_(Foam::constant::mathematical::pi),
    anglemin_(-1*Foam::constant::mathematical::pi),
    omegamax_(this->anglemax_),
    omegamin_(this->anglemin_)
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
    controlDelay_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("controlDelay",1.0)),
    controlTarget_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault<int>("controlTarget",2)),
    forcesDict_(PIDcontrolDict_.subDict("PIDcontroller").subDict("controlledVarData")),
    setPoint_(forcesDict_.getOrDefault<scalar>("setPoint",0.35)),
    direction_(forcesDict_.getOrDefault<label>("direction",3)),
    mycoefs(6,Zero),
    myForce(Zero),
    anglemax_(Foam::constant::mathematical::pi),
    anglemin_(-1*Foam::constant::mathematical::pi),
    omegamax_(this->anglemax_),
    omegamin_(this->anglemin_)
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
    // esto no lo lee el hijoputa
    //Info<<"AQUI TAMPOCO?"<<endl;
    //PIDcontrolDict_.
    Info<<"P: "<<P_<<endl;
    Info<<"I: "<<I_<<endl;
    Info<<"D: "<<D_<<endl;
    Info<<"setPoint: "<<setPoint_<<endl;
    Info<<"Control start delay: "<<controlDelay_<<endl;
    Info<<"Full PIDcontrolDict file:"<<endl;
    Info<<"---------------------------------------------"<<endl;
    Info<<PIDcontrolDict_<<endl;
    if (PIDcontrolDict_.subDict("PIDcontroller").isDict("actuatorModelCoeffs"))
    {
        //Info<< "GOT HERE"<<endl;
        dictionary actuatorModel(PIDcontrolDict_.subDict("PIDcontroller").subDict("actuatorModelCoeffs"));
        anglemax_ = actuatorModel.getOrDefault<scalar>("angleMax",anglemax_);
        anglemin_ = actuatorModel.getOrDefault<scalar>("angleMin",anglemin_);
        omegamax_ = actuatorModel.getOrDefault<scalar>("omegaMax",omegamax_);
        omegamin_ = actuatorModel.getOrDefault<scalar>("omegaMin",omegamin_);
    }
    
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
    controlDelay_(ptf.controlDelay_),
    controlTarget_(ptf.controlTarget_),
    forcesDict_(ptf.forcesDict_),
    setPoint_(ptf.setPoint_),
    direction_(ptf.direction_),
    mycoefs(6,Zero),
    myForce(Zero),
    anglemax_(ptf.anglemax_),
    anglemin_(ptf.anglemin_),
    omegamax_(ptf.omegamax_),
    omegamin_(ptf.omegamin_)
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
    controlDelay_(ptf.controlDelay_),
    controlTarget_(ptf.controlTarget_),
    forcesDict_(ptf.forcesDict_),
    setPoint_(ptf.setPoint_),
    direction_(ptf.direction_),
    mycoefs(6,Zero),
    myForce(Zero),
    anglemax_(ptf.anglemax_),
    anglemin_(ptf.anglemin_),
    omegamax_(ptf.omegamax_),
    omegamin_(ptf.omegamin_)
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
    // PENDING TO DECLARE CONTROLDELAY VAR
    if (t.value() < controlDelay_)
    {
        angle = angle0_ + omega_*t.value();
        vector axisHat = axis_/mag(axis_);
        vectorField p0Rel(p0_ - origin_);
        // report calculated angles
        Info<<"Control not activated yet"<<endl;
        report(Info);
        vectorField::operator=
        (
            p0Rel*(cos(angle) - 1)
        + (axisHat ^ p0Rel*sin(angle))
        + (axisHat & p0Rel)*(1 - cos(angle))*axisHat
        );

        fixedValuePointPatchField<vector>::updateCoeffs();

        return;
    }
    
    scalar dt(t.deltaTValue());

    if (timeIndex_!= t.timeIndex())
    {
        timeIndex_ = t.timeIndex();
        oldOmega_= omega_;
        oldError_ = error_;
        oldErrorIntegral_ = errorIntegral_;
        oldangle_ = angle;
    }

    FIFOStack<scalar> windowTimes; // work in progress

    //functionObjects::lforces f("forces",t,forcesDict_,true);
    functionObjects::lforceCoeffs fc("forceCoeffs",t,forcesDict_,true);

    switch (controlTarget_) {
        case 0:
            // not implemented yet
            
            //functionObjects::forces f("forces", mesh, forcesDict_);
            fc.calcForcesMoment();

            myForce = fc.forceEff(); // THIS IS STILL PENDING
            error_ = setPoint_ - myForce.x(); // this should choose x,y,z according to direction

            errorIntegral_ = oldErrorIntegral_ + I_*0.5*(error_ + oldError_)*dt;
            errorDifferential_ = (error_ - oldError_)/dt;
            //const scalar errorDifferential = oldError_ - error_;
            angle = P_*error_ + errorIntegral_ + D_*errorDifferential_;
            omega_ = (angle - oldangle_)/dt;

            Info<< "Controlled variable (Lift) : "<< myForce[2]<<endl;
            //Info<< "dir1 force is: "<< myForce.x() <<endl;
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
            //Info<<"forceCoeffs"<<endl;
            //fc.calcForcesMoment();
            fc.calcrot();
            mycoefs = fc.coefList;
            //Info<<"Cd"<< mycoefs[0]<<endl;
            //Info<<"Cs"<< mycoefs[1]<<endl;
            Info<<"Controlled Variable (Cl): "<< mycoefs[2]<<endl;
            //Info<<"CmRoll"<< mycoefs[3]<<endl;
            //Info<<"CmPitch"<< mycoefs[4]<<endl;
            //Info<<"CmYaw"<< mycoefs[5]<<endl;

            error_ = setPoint_ - mycoefs[2]; // pending option to choose other coeffs than Cl
            errorIntegral_ = oldErrorIntegral_ + I_*0.5*(error_ + oldError_)*dt;
            errorDifferential_ = (error_ - oldError_)/dt;
            angle = P_*error_ + errorIntegral_ + D_*errorDifferential_;
            omega_ = (angle - oldangle_)/dt; 

            

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

            error_ = setPoint_ - angle;

            errorIntegral_ = oldErrorIntegral_ + I_*0.5*(error_ + oldError_)*dt;
            errorDifferential_ = (error_ - oldError_)/dt;
            angle = P_*error_ + errorIntegral_ + D_*errorDifferential_;
            omega_ = (angle - oldangle_)/dt;

            break;

        default:

            Info<< "Unknown controlled variable"<<endl
            <<"Allowed variables are:"<< endl<< "(forces forceCoeffs angle)"
            <<endl<<exit(FatalError);

        break;

        // EXTERNAL CONTROL WILL BE DONE IN ANOTHER pointPatchFields option
    }

    rawomega_ = omega_;
    rawangle_ = angle;

    // SATURATOR. LIMITS MAXIMUM ANGLE & OMEGA

    // 1. Check omega limits and modify itself and angle accordingly
    // REQUIRES CORRECT OMEGA CALCULATION IN ABOVE CODE LINES!!
    // you could declare a state vector that contains angle and omega so that a function can be called
    // you should learn how to use the tensor and rotation classes to keep track of these things??

    if (pos(omega_) && omega_>omegamax_)
    {
        Info<<"Limiting omega "<<endl;//from "<< omega_<< " to "<< omegamax_<< " (rad/s)"<<endl;
        //Info<<"Control signal was "<< angle;
        omega_ = omegamax_;
        angle = oldangle_ + omega_*dt;
        //Info<<"(rad). Corrected signal is "<<angle<<"(rad)"<<endl;
    }
    else if (neg(omega_) && omegamin_>omega_)
    {
        Info<<"Limiting omega "<<endl;//from "<< omega_<< " to "<< omegamax_<< " (rad/s)"<<endl;
        //Info<<"Control signal was "<< angle;
        omega_ = omegamin_;
        angle = oldangle_ + omega_*dt;
        //Info<<"(rad). Corrected signal is "<<angle<<"(rad)"<<endl;
    }

    if (pos(angle) && angle>anglemax_)
    {
        Info<<"Limiting angle"<<endl;// from "<<angle<<" to "<<anglemax_<<" (rad)"<<endl;
        //Info<<"Mesh omega was "<<omega_;
        omega_ = (anglemax_ - oldangle_)/dt;
        angle = anglemax_;
        //Info<<" (rad/s). Corrected mesh omega is "<<omega_<<"(rad/s)"<<endl;
    }
    else if (neg(angle) && angle<anglemin_)
    {
        Info<<"Limiting angle"<<endl;// from "<<angle<<" to "<<anglemin_<<" (rad)"<<endl;
        //Info<<"Mesh omega was "<<omega_;
        omega_ = (anglemin_ - oldangle_)/dt;
        angle = anglemin_;
        //Info<<" (rad/s). Corrected mesh omega is "<<omega_<<"(rad/s)"<<endl;
    }
    //scalar signo(sign(-error_));
    //Info<< "signo es: "<<signo<<endl;
    // MOVEMENT BASED ON CALCULATED OMEGA
    //angle = angle0_ + angle ;//+ omega_*dt;
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
    os.writeEntry("omegamin", omegamin_);
    p0_.writeEntry("p0", os);
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
            IOobject::MUST_READ,   // PENDING TO WORK ON A VARIABLE OR PRESCRIBED SETPOINT 
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
    os.writeEntry("Raw omega", rawomega_);
    os.writeEntry("Raw angle",rawangle_);
    os.writeEntry("Saturated omega", omega_);
    os.writeEntry("Saturated angle", angle);
    os.writeEntry("error", error_);
    os.writeEntry("errorIntegral", errorIntegral_);
    os.writeEntry("errorDifferential",errorDifferential_);

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
