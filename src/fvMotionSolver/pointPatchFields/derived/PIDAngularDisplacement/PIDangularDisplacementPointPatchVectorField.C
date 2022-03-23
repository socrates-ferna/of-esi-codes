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
#include "FIFOStack.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Enum<PIDangularDisplacementPointPatchVectorField::filterType>
PIDangularDisplacementPointPatchVectorField::filterTypeNames
({
    { filterType::FIRST, "firstOrderSimple" },
    { filterType::SECOND, "secondOrderSimple" },
    { filterType::BUTTER, "butterworth"},
    { filterType::AVERAGE, "average"}
});

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
    angle(this->angle0_),
    PIDcontrolDict_(this->readControl()), 
    P_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("P",0.5)),
    I_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("I",0.5)),
    D_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("D",1.0)),
    T_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("T",1e-4)),
    //PIDfilter_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault<bool>("COfilter",false)),
    errorIntegral_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("errorIntegral",0.0)),
    errorDifferential_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("errorDifferential",0.0)),
    error_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("error",0.0)),
    controlDelay_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault("controlDelay",1.0)),
    //applyFilter_(false),
    controlTarget_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault<int>("controlTarget",2)),
    forcesDict_(PIDcontrolDict_.subDict("PIDcontroller").subDict("controlledVarData")),
    //fc_d("forceCoeffs",db().time(),forcesDict_),
    setPoint_(forcesDict_.getOrDefault<scalar>("setPoint",0.35)),
    direction_(forcesDict_.getOrDefault<scalar>("direction",3)),
    mycoefs(6,Zero),
    myForce(Zero),
    //saturate_(false),
    anglemax_(Foam::constant::mathematical::pi),
    anglemin_(-1*Foam::constant::mathematical::pi),
    omegamax_(this->anglemax_),
    omegamin_(this->anglemin_),
    controlString_(this->controlString(forcesDict_,controlTarget_)),
    antiWindup_(PIDcontrolDict_.subDict("PIDcontroller").get<bool>("antiwindup")),
    integralReset_(PIDcontrolDict_.subDict("PIDcontroller").get<bool>("integralReset")),

    //meanValue_(forcesDict_.getOrDefault<scalar>("initialValue",0.0)),
    //PVfilterDict_(PIDcontrolDict_.subDict("PIDcontroller").subDict("PVfilter")),
    //filterType_(0),
    //window_(PVfilterDict_.getOrDefault<scalar>("window",0.01)),
    t0_(Zero)/*,
    testForces("forceCoeffs",mesh.time(),forcesDict_,true)*/
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
    p0_(p.localPoints()),
    angle(this->angle0_),
    PIDcontrolDict_(this->readControl()),
    P_(PIDcontrolDict_.subDict("PIDcontroller").get<scalar>("P")),
    I_(PIDcontrolDict_.subDict("PIDcontroller").get<scalar>("I")),
    D_(PIDcontrolDict_.subDict("PIDcontroller").get<scalar>("D")),
    T_(PIDcontrolDict_.subDict("PIDcontroller").get<scalar>("T")),
    //PIDfilter_(PIDcontrolDict_.subDict("PIDcontroller").getOrDefault<bool>("COfilter",false)),
    errorIntegral_(PIDcontrolDict_.subDict("PIDcontroller").get<scalar>("errorIntegral")),
    errorDifferential_(PIDcontrolDict_.subDict("PIDcontroller").get<scalar>("errorDifferential")),
    error_(PIDcontrolDict_.subDict("PIDcontroller").get<scalar>("error")),
    controlDelay_(PIDcontrolDict_.subDict("PIDcontroller").get<scalar>("controlDelay")),
    //applyFilter_(false),
    controlTarget_(PIDcontrolDict_.subDict("PIDcontroller").get<int>("controlTarget")), //DO IT WITH ENUM?
    forcesDict_(PIDcontrolDict_.subDict("PIDcontroller").subDict("controlledVarData")),
    //fc_d("forceCoeffs",db().time(),forcesDict_),
    setPoint_(forcesDict_.get<scalar>("setPoint")),
    direction_(forcesDict_.get<label>("direction")),
    mycoefs(6,Zero),
    myForce(Zero),
    //saturate_(false),
    anglemax_(Foam::constant::mathematical::pi),
    anglemin_(-1*Foam::constant::mathematical::pi),
    omegamax_(this->anglemax_),
    omegamin_(this->anglemin_),
    controlString_(this->controlString(forcesDict_,controlTarget_)),
    meanValue_(forcesDict_.getOrDefault<scalar>("initialValue",0.0)),
    antiWindup_(PIDcontrolDict_.subDict("PIDcontroller").get<bool>("antiwindup")),
    integralReset_(PIDcontrolDict_.subDict("PIDcontroller").get<bool>("integralReset")),
    //REMOVE.CONFLICTS WITH PV_
    //PVfilterDict_(PIDcontrolDict_.subDict("PIDcontroller").subDict("PVfilter")), //REMOVE AND DO IN INITALISATION
    t0_(dict.getOrDefault<scalar>("t0",0.0))/*,
    testForces("forceCoeffs",mesh.time(),forcesDict_,true)*/
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
    
    if (dict.found("p0"))
    {
        p0_ = vectorField("p0", dict , p.size());
        //angle0_ = dict.get<scalar>("angle");
        //Info<< "RESETTING ANGLE: "<< angle0_<<endl;
    }
    /*else
    {
        p0_ = p.localPoints();
    }
    */
    p0_ = p.localPoints();
    Info<<"P: "<<P_<<endl; // PONLO TODO EN UNA LINEA
    Info<<"I: "<<I_<<endl;
    Info<<"D: "<<D_<<endl;
    Info<<"setPoint: "<<setPoint_<<endl;
    Info<<"Control start delay: "<<controlDelay_<<endl;
    Info<<"Full PIDcontrolDict file:"<<endl;
    Info<<"---------------------------------------------"<<endl;
    Info<<PIDcontrolDict_<<endl;

    // Actuator model setup. Saturation + roll-off
    if (PIDcontrolDict_.subDict("PIDcontroller").isDict("actuatorModelCoeffs"))
    {
        actuatorModel = PIDcontrolDict_.subDict("PIDcontroller").subDict("actuatorModelCoeffs");
        saturate_ = actuatorModel.getOrDefault<bool>("saturate",false);
        
        if(saturate_)
        {
        anglemax_ = actuatorModel.getOrDefault<scalar>("angleMax",anglemax_);
        anglemin_ = actuatorModel.getOrDefault<scalar>("angleMin",anglemin_);
        omegamax_ = actuatorModel.getOrDefault<scalar>("omegaMax",omegamax_);
        omegamin_ = actuatorModel.getOrDefault<scalar>("omegaMin",omegamin_);
        }
        else {Warning<<" NO SATURATION APPLIED TO CONTROL SIGNAL"<<endl;}

        COrolloff_ = actuatorModel.getOrDefault<bool>("rolloff",false);
        if (COrolloff_){
            // Both instructions should raise fatal error if keywords not correctly specified
            COrolloffType_ = filterTypeNames.get("rolloffType",actuatorModel);
            // FALTA UN IF POR SI ALGÚN LOCO PONE AQUÍ UN BUTTERWORTH
            int n = static_cast<int>(COrolloffType_);
            Info<<"Control Signal roll-off. Filter Type: "<<n<<endl;
            selectFilterAndSetup(COnum_,COdenom_,COrolloffType_,actuatorModel, dict);
            // PENDING TO CALL FILTER SETUP
        }
        else {Info<<"No roll-off applied to Control Signal"<<endl;}
    }
    else {Warning<<" NO ACTUATOR MODEL OR ROLL-OFF FILTER SPECIFIED. RISK OF HIGH MESH DEFORMATION OR dt->0 to control MeshCo if pimpleDyMStopFoam is used"<<endl;}

    // Feedback filtering setup
    if (PIDcontrolDict_.subDict("PIDcontroller").isDict("PVfilter"))
    {
        PVfilterDict_ = PIDcontrolDict_.subDict("PIDcontroller").subDict("PVfilter");
        applyFilter_ = PVfilterDict_.getOrDefault<bool>("applyFilter",false);
        if(applyFilter_)
        {
            PVfilterType_ = filterTypeNames.get("type",PVfilterDict_);

            int n = static_cast<int>(PVfilterType_);
            Info<<"Feedback Signal filtering. Filter Type: "<<n<<endl;

            selectFilterAndSetup(PVnum_, PVdenom_, PVfilterType_, PVfilterDict_, dict);
        }
        else {Warning<<"NO FILTER APPLIED TO PROCESS VARIABLE SIGNAL. RISK OF NOISY FEEDBACK"<<endl;}
    }
    else {Warning<<"NO PROCESS VARIABLE FILTER DICTIONARY AVAILABLE. RISK OF NOISY FEEDBACK"<<endl;}

    PVinbuffer_ = {0,0,0}; COinbuffer_ = {0,0,0};
    PVoutbuffer_ = {0,0,0}; COoutbuffer_ = {0,0,0};

    //bool test = PIDcontrolDict_.subDict("PIDcontroller").getBool("test");

    if (PIDcontrolDict_.subDict("PIDcontroller").getOrDefault<bool>("test",false))
    {
        // Perform test for step input to digital filters
        scalar stepInput = 1.0;
        PV_ = stepInput;
        angle = stepInput;
        scalar filteredPV;
        scalar filteredAngle;
        scalar simTime = 0.0;
        scalar simFinalTime = 0.1;
        scalar localdt = T_;
        int i=0;
        while(simTime<simFinalTime){
    //const vector& num, const vector& denom, vector& inbuf, vector& outbuf, const filterType filterName, scalar& dt //revisa q los const valgan de algo así en el libro
            //PVfilter call
            Info<<i<<"        "<<endl;
            if (applyFilter_){
                filteredPV = digitalFilter(PV_,PVnum_,PVdenom_,PVinbuffer_,PVoutbuffer_,PVfilterType_,localdt);
                Info<<"PV:"<<filteredPV<<endl;

            }
            //COfilter call
            if (COrolloff_){
                filteredAngle = digitalFilter(angle,COnum_,COdenom_,COinbuffer_,COoutbuffer_,COrolloffType_,localdt);
                Info<<"CO:"<<filteredAngle<<endl;
            }
            simTime+=localdt;
            i++;
        }
        Info<<"simTime:"<<simTime<<endl;
        Info<<"PV: "<<PV_<<endl;
        Info<<"CO: "<<angle<<endl;
        Info<<"filtered PV: "<<filteredPV<<endl;
        Info<<"filtered angle: "<<filteredAngle<<endl;
        Info<<"Iterations: "<<i<<endl;
        Info<<exit(FatalError);
    }
    
    //Info<<exit(FatalError);
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
    angle(ptf.angle0_),
    PIDcontrolDict_(ptf.PIDcontrolDict_),
    P_(ptf.P_),
    I_(ptf.I_),
    D_(ptf.D_),
    T_(ptf.T_),
    //PIDfilter_(ptf.PIDfilter_),
    errorIntegral_(ptf.errorIntegral_),
    errorDifferential_(ptf.errorDifferential_),
    error_(ptf.error_),
    controlDelay_(ptf.controlDelay_),
    applyFilter_(ptf.applyFilter_),
    controlTarget_(ptf.controlTarget_),
    forcesDict_(ptf.forcesDict_),
    //fc_d(ptf.fc_d),
    setPoint_(ptf.setPoint_),
    direction_(ptf.direction_),
    mycoefs(6,Zero),
    myForce(Zero),
    saturate_(ptf.saturate_),
    anglemax_(ptf.anglemax_),
    anglemin_(ptf.anglemin_),
    omegamax_(ptf.omegamax_),
    omegamin_(ptf.omegamin_),
    controlString_(ptf.controlString_),
    antiWindup_(ptf.antiWindup_),
    integralReset_(ptf.integralReset_),
    //meanValue_(ptf.meanValue_),
    //PVfilterDict_(ptf.PVfilterDict_),
    //filterType_(ptf.filterType_),
    //window_(ptf.window_),
    t0_(ptf.t0_)/*,
    testForces("forceCoeffs",mesh.time(),forcesDict_,true)*/
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
    angle(ptf.angle0_),
    PIDcontrolDict_(ptf.PIDcontrolDict_),
    P_(ptf.P_),
    I_(ptf.I_),
    D_(ptf.D_),
    T_(ptf.T_),
    //PIDfilter_(ptf.PIDfilter_),
    errorIntegral_(ptf.errorIntegral_),
    errorDifferential_(ptf.errorDifferential_),
    error_(ptf.error_),
    controlDelay_(ptf.controlDelay_),
    applyFilter_(ptf.applyFilter_),
    controlTarget_(ptf.controlTarget_),
    forcesDict_(ptf.forcesDict_),
    //fc_d(ptf.fc_d),
    setPoint_(ptf.setPoint_),
    direction_(ptf.direction_),
    mycoefs(6,Zero),
    myForce(Zero),
    saturate_(ptf.saturate_),
    anglemax_(ptf.anglemax_),
    anglemin_(ptf.anglemin_),
    omegamax_(ptf.omegamax_),
    omegamin_(ptf.omegamin_),
    controlString_(ptf.controlString_),
    antiWindup_(ptf.antiWindup_),
    integralReset_(ptf.integralReset_),
    //meanValue_(ptf.meanValue_),
    //PVfilterDict_(ptf.PVfilterDict_),
    //filterType_(ptf.filterType_),
    //window_(ptf.window_),
    t0_(ptf.t0_)/*,
    testForces("forceCoeffs",mesh.time(),forcesDict_,true)*/
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
    if (t.value() < controlDelay_)
    {
        if (timeIndex_!=t.timeIndex())
        {
            timeIndex_=t.timeIndex();
        }
        angle = angle0_ + omega_*t.value(); // you should change it back to oscillatory
        scalar passedAngle = angle - angle0_;
        vector axisHat = axis_/mag(axis_);
        vectorField p0Rel(p0_ - origin_);
        // report calculated angles
        Info<<"Control not activated yet"<<endl;
        Info<<"Time index: "<< t.timeIndex()<<endl;
        report(Info);
        vectorField::operator=
        (
            p0Rel*(cos(passedAngle) - 1)
        + (axisHat ^ p0Rel*sin(passedAngle))
        + (axisHat & p0Rel)*(1 - cos(passedAngle))*axisHat
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
        oldangle_ = angle ;
    }


    functionObjects::lforceCoeffs fc("forceCoeffs",t,forcesDict_,true);

    switch (controlTarget_) {
        case 0:

            fc.calcForcesMoment();
            myForce = fc.forceEff(); // forceEff from lforces Class RETURNS ROTATED FORCE
            PV_=myForce[direction_];
            //rollavg(myForce[direction_],PVfilterDict_,dt); // rolling Average with FIFOStack class
            if(applyFilter_){
                PV_ = digitalFilter(PV_,PVnum_, PVdenom_,PVinbuffer_,PVoutbuffer_,PVfilterType_,dt);
                Info<<"Unfiltered PV from function is: "<<myForce[direction_]<<" ("<<controlString_<<")"<<endl;
                Info<<"Filtered PV from function is: "<<PV_<<" ("<<controlString_<<")"<<endl;
            }
            //Calculate unfiltered Control Output
            rawCO(dt,PV_);
            
            Info<< "Controlled variable"<<controlString_<<" : "<< myForce[direction_]<<endl;

        break;

        case 1: 

            fc.calcrot();
            mycoefs = fc.coefList;
            PV_=mycoefs[direction_];
            
            if (applyFilter_) {
                PV_ = digitalFilter(PV_,PVnum_, PVdenom_,PVinbuffer_,PVoutbuffer_,PVfilterType_,dt);
                Info<<"Unfiltered PV from function is: "<<mycoefs[direction_]<<" ("<<controlString_<<")"<<endl;
                Info<<"Filtered PV from function is: "<<PV_<<" ("<<controlString_<<")"<<endl;
            }
            rawCO(dt,PV_);

            break;

        case 2:

            rawCO(dt,angle);
            break;

        default:
            Info<< "Unknown controlled variable"<<endl
            <<"Allowed variables are:"<< endl<< "(forces forceCoeffs angle)"
            <<endl<<exit(FatalError);

        break;
    }

    if(COrolloff_) {
        angle = digitalFilter(angle,COnum_, COdenom_, COinbuffer_,COoutbuffer_,COrolloffType_,dt);
    }

    omega_ = (angle - oldangle_)/dt;
    rawomega_ = omega_;
    rawangle_ = angle;

    if(saturate_) {
        saturator(dt);
    }
    satDiff_ = angle - rawangle_;
    

    vector axisHat = axis_/mag(axis_);
    vectorField p0Rel(p0_ - origin_);
    // report calculated angles
    report(Info);
    scalar passedAngle = angle - angle0_;
    vectorField::operator=
    (
        p0Rel*(cos(passedAngle) - 1)
      + (axisHat ^ p0Rel*sin(passedAngle))
      + (axisHat & p0Rel)*(1 - cos(passedAngle))*axisHat
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
    //Info<< "db access"<< db().names() << endl;

}


void PIDangularDisplacementPointPatchVectorField::rawCO(const scalar dt, const scalar PVsignal)
{
    scalar olderror = error_;
    error_ = setPoint_ - PVsignal;
    errorIntegral_ = oldErrorIntegral_ + I_*0.5*(error_ + oldError_)*dt;
    if (antiWindup_){
        errorIntegral_ = errorIntegral_ + I_/P_*satDiff_*dt;
        Info<<"SatDiff= "<<satDiff_<<" Integral action= "<< I_*0.5*(error_ + oldError_)*dt <<" Antiwindup action= "<< I_/P_*satDiff_*dt<<endl;
    }
    errorDifferential_ = (error_ - oldError_)/dt;
    angle = P_*error_ + errorIntegral_ + D_*errorDifferential_;
    if (integralReset_){
        if (neg(olderror*error_)){
            Info<<"Error changed sign, Integral action reset to 0"<<endl;
            errorIntegral_ = 0;
        }
    }
    omega_ = (angle - oldangle_)/dt;
    //Info<<"rawCO called, Output: "<<angle<<" raw omega: "<<omega_<<endl;
}

void PIDangularDisplacementPointPatchVectorField::saturator(const scalar dt)
{
        // SATURATOR. LIMITS MAXIMUM ANGLE & OMEGA

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
    os.writeEntry("t0",t0_);
    os.writeEntry("P",P_);
    os.writeEntry("I",I_);
    os.writeEntry("D",D_);
    if (PVfilterDict_.found("windowValues")) // OJO AQUÍ QUE PVfilterDict no tiene por qué ser
    {
        os.writeEntry("windowValues",PVfilterDict_.get<scalarList>("windowValues"));
    }
    if (PVfilterDict_.found("windowTimes"))
    {
        os.writeEntry("windowTimes",PVfilterDict_.get<scalarList>("windowTimes"));
    }
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
            IOobject::MUST_READ_IF_MODIFIED,   // PENDING TO WORK ON A VARIABLE OR PRESCRIBED SETPOINT 
            IOobject::AUTO_WRITE
        )
    );
    
    Info<< "Read controller setup from PIDcontrolDict:"<< endl<< readdict.dictName() << endl;
    return readdict;
}

// Determine control variable string to output during execution
word PIDangularDisplacementPointPatchVectorField::controlString(const dictionary& dict,label controlTarget)
{
        wordList forceStrings = {"Drag Force","Side Force", "Lift Force", "Roll Moment", "Pitch Moment", "Yaw Moment"};
        wordList forceCoeffStrings = {"Cd", "Cs", "Cl", "Cr", "Cp", "Cy"};
        
        //Info<<forceStrings<<endl;
        //Info<<forceCoeffStrings<<endl;
        
        word controlString = "Cl";

        label direc = dict.get<int>("direction");

        if (controlTarget == 0)
        {
            controlString = forceStrings[direc];
        }
        else if (controlTarget == 1)
        {
            controlString = forceCoeffStrings[direc];
        }
        else
        {
            controlString = "angle";
        }

        return controlString;
}

// CURRENTLY ONLY AVAILABLE FOR SCALARS. 
//IDEAL WOULD BE TO IMITATE TEMPLATED FUNCTION AS IN 
//averageConditionTemplates.C
void PIDangularDisplacementPointPatchVectorField::rollavg
(
    scalar currentValue,
    dictionary& dict,
    const scalar& dt
)
{
    FIFOStack<scalar> windowTimes;
    FIFOStack<scalar> windowValues;

    dict.readIfPresent("windowTimes",windowTimes);
    dict.readIfPresent("windowValues",windowValues);

    for (scalar& dti : windowTimes) //if windowTimes empty will be zero?
    {
        dti += dt + t0_; //increase value of all elements by dt
    }

    bool removeValue = true;
    // check if values in stack have reached the window
    while (removeValue && windowTimes.size())
    {
        removeValue = (windowTimes.first() - t0_) > window_;

        if (removeValue)
        {
            windowTimes.pop();
            windowValues.pop();
        }
    }

    // Add the current value. Up to here nothing must have happenned if it's empty
    windowTimes.push(dt);
    windowValues.push(currentValue);
    
    // Define iterators Calculate the average
    auto timeIter = windowTimes.cbegin();
    auto valueIter = windowValues.cbegin();

    scalar meanValue = 0; // pTraits<Type>::zero in averageConditionTemplates.C
    scalar valueOld(0);

    for
    (
        label i = 0;
        timeIter.good();
        ++i,++timeIter,++valueIter
    )
    {
        const scalar& value = valueIter();
        const scalar dt = timeIter();

        meanValue += dt*value;

        if(i){ meanValue -= dt*valueOld;}

        valueOld = value;
    }

    meanValue /= windowTimes.first();

    // Store the state informtion for the next call
    dict.set("windowTimes",windowTimes);
    dict.set("windowValues",windowValues);

    PV_ = meanValue;
    scalar delta = meanValue - currentValue;
    
    Info<<"meanValue is: "<<PV_<<"  Delta is: "<<delta<<endl;
    /*
    Info<<"FIFO Window Times "<<windowTimes<<endl;
    Info<<"FIFO Window Values"<<windowValues<<endl;
    */
}

scalar PIDangularDisplacementPointPatchVectorField::digitalFilter
(
    scalar inSignal, const vector& num, const vector& denom, vector& inbuf, vector& outbuf, const filterType filterName, scalar& dt //revisa q los const valgan de algo así en el libro
)
{
    // PUEDO MONTAR DOS BUFFERS E IRLOS ACTUALIZANDO EN CADA ITERACIÓN
    // DE ESA FORMA SOLO ME HACE FALTA ESTA FUNCIÓN PARA LAS FÓRMULAS EN DIFERENCIAS
    // inbuf[0]=x(k), inbuf[1] = x(k-1), inbuf[2] = x(k-2)
    // outbuf[0]=y(k-1), outbuf[1] = y(k-2), outbuf[2] = y(k-3) // We enter the loop with outbuf from previous it (outbut[2] not necessary)
    scalar outSignal = 0.0;

    inbuf[2]=inbuf[1];inbuf[1]=inbuf[0];inbuf[0]=inSignal;
    Info<<"Inbuf:"<<inbuf;
    Info<<"     Outbuf"<<outbuf;

    switch (filterName)
    {
    case filterType::BUTTER :
        outSignal = num[0]*inbuf[2] + num[1]*inbuf[1] + num[2]*inbuf[0] 
        - denom[0]*outbuf[1] - denom[1]*outbuf[0]; //- denom[2]*outbuf[0];
        break;
    case filterType::FIRST :
        outSignal = num[0]*inbuf[1] + num[1]*inbuf[0] - denom[0]*outbuf[0];// - denom[1]*outbuf[0];
        break;
    case filterType::SECOND :
        outSignal = num[0]*inbuf[2] + num[1]*inbuf[1] + num[2]*inbuf[0] 
        - denom[0]*outbuf[1] - denom[1]*outbuf[0];// - denom[2]*outbuf[0];
        break;
    case filterType::AVERAGE :
        rollavg(outSignal,PVfilterDict_,dt);
        break;
    default:
        Info<<"Unknown filter type"<<exit(FatalError)<<endl;
        break;
    }
    outbuf[2]=outbuf[1];outbuf[1]=outbuf[0];outbuf[0]=outSignal;
    Info<<"    After filter:"<<outbuf<<endl;
    // slip
    return outSignal;
}

PIDangularDisplacementPointPatchVectorField::tf PIDangularDisplacementPointPatchVectorField::zplane_tf
(
    vector numerator, vector denominator // debería pasarlas por referencia y olvidarme del
)
{
    // Num and denom are passed as continuous transfer function coefs bn,an
    tf digital_tf;
    // prewarping should be done before this function directly to wc
    scalar K = 2/T_;
    scalar Ksq = K*K;
    
    
    float bz0, bz1, bz2, az0, az1, az2;
    if ((denominator[2] == 0) && (numerator[2] == 0)){
        Info<<"First order Filter"<<endl;
        //         b0*s + b1
        // G(s) = -----------
        //         a0*s + a1

        //         bz0*z^-1 + bz1
        // G(z) = ---------------------------
        //         az0*z^-1 + az1
        // denom = a0*K + a1
        // bz0 = (-b0*K + b1) / denom
        // bz1 = (b0*K + b1) / denom
        // az0 = (-a0*K + a1) / denom
        // az1 = 1
        scalar denom = denominator[0]*K + denominator[1];
        bz2 = 0; az2 = 0;
        bz0 = (-numerator[0]*K + numerator[1]) / denom;
        bz1 = (numerator[0]*K + numerator[1]) / denom;
        az0 = (-denominator[0]*K + denominator[1]) / denom;
        az1 = 1.0;
        digital_tf.denominator = {az0,az1,az2};
        digital_tf.numerator = {bz0,bz1,bz2};

    }
    else {
        //         b0*s^2 + b1*s + b2
        // G(s) = -------------------
        //         a0*s^2 + a1*s + a2

        //         bz0*z^-2 + bz1*z^-1 + bz2
        // G(z) = ---------------------------
        //         az0*z^-2 + az1*z^-1 + az2

        // denom = a2/K^2 + a1/K + a0
        // bz0 = (b2/K^2 - b1/K + b0) / denom
        // bz1 = (2*b2/K^2 - 2*b0) /denom
        // bz2 = (b2/K^2 + b1/K + b0) /denom
        // same for a_n coefs
        // **az2 is 1 because G(z) is normalised by az2 for convenience in the difference equation (therefore the division by denom)
        Info<<"Second order Filter"<<endl;
        scalar denom = denominator[2]/Ksq + denominator[1]/K + denominator[0];
        bz0 = (numerator[2]/Ksq - numerator[1]/K + numerator[0]) / denom;
        bz1 = (2*numerator[2]/Ksq - 2*numerator[0]) / denom;
        bz2 = (numerator[2]/Ksq + numerator[1]/K + numerator[0]) / denom;
        az0 = (denominator[2]/Ksq - denominator[1]/K + denominator[0]) / denom;
        az1 = (2*denominator[2]/Ksq - 2*denominator[0]) / denom;
        az2 = 1.0;

        digital_tf.numerator = {bz0,bz1,bz2};
        digital_tf.denominator = {az0,az1,az2};

    }
    
    Info<<"Transfer function G(s)="<<numerator<<" / "<<denominator<<endl;
    Info<<"Transfer function G(z)="<<digital_tf.numerator<<" / "<<digital_tf.denominator<<endl;
    
    return digital_tf;
}

void PIDangularDisplacementPointPatchVectorField::selectFilterAndSetup
(
    vector& tfnum, vector& tfdenom, const filterType filterName, 
    dictionary& filterDict, const dictionary& pdict
)
{
    tf tftf;
    scalar wc;
    scalar Tf;
    switch(filterName){
        case filterType::BUTTER :
            wc = filterDict.get<scalar>("wc");
            wc = 2/T_*tan(wc*T_/2); // pre-warping. wc ~ wcwarped if |wc|<<pi/T
            tfnum = {0,0,wc*wc};
            tfdenom = {1,wc*sqrt(2.0),wc*wc};
            tftf = zplane_tf(tfnum,tfdenom);
            tfnum = tftf.numerator;
            tfdenom = tftf.denominator;
            break;
        case filterType::FIRST :
            Tf = filterDict.get<scalar>("Tf");
            tfnum = {0,1,0};
            tfdenom = {Tf,1,0};
            tftf = zplane_tf(tfnum,tfdenom);
            tfnum = tftf.numerator;
            tfdenom = tftf.denominator;
            break;
        case filterType::SECOND :
            Tf = filterDict.get<scalar>("Tf");
            tfnum = {0,0,1};
            tfdenom = {Tf*Tf,2*Tf,1};
            tftf = zplane_tf(tfnum,tfdenom);
            tfnum = tftf.numerator;
            tfdenom = tftf.denominator;
            break;
        case filterType::AVERAGE :
            window_ = filterDict.get<scalar>("window");
            PV_ = forcesDict_.getOrDefault<scalar>("initialValue",0.0);
            if (pdict.found("windowValues")) // aquí puedes usar el dict tal como escribe la funcion report!!!
            {
                scalarList catchWindowValues(pdict.get<scalarList>("windowValues"));
                filterDict.set("windowValues",catchWindowValues);
                Info<<"READ WINDOWVALUES"<<endl;
            }
            if (pdict.found("windowTimes"))
            {
                scalarList catchWindowTimes(pdict.get<scalarList>("windowTimes"));
                filterDict.set("windowValues",catchWindowTimes);
                Info<<"READ WINDOWTIMES"<<endl;
            }
            break;
        default :
            Info<<"Unknown Filter Type"<<exit(FatalError);
            break;                   
    }
    return;
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
    os.writeEntry("Angle seen by OpenFOAM",angle-angle0_);
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
