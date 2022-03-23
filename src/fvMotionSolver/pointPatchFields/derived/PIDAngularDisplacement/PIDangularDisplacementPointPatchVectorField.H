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
#include "lforces.H"
#include "lforceCoeffs.H"

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
    public:

    // filterType class declaration ahead of
    enum class filterType{
            FIRST,
            SECOND,
            BUTTER,
            AVERAGE
        };
    static const Enum<filterType> filterTypeNames;

    // Transfer function. Serves as container for denom and num vectors result from zplane_tf()
    struct tf{vector numerator; vector denominator;};


    private: //Or protected?
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
        // Control Output Smoothing through roll-off filter, e.g. G(s) = 1 / (Tf*s + 1)^2
        // We use it because optiPID script from Octave control pkg for PID tuning assumes roll-off function
        // It is similar to an actuator connected to the raw control signal
        bool COrolloff_=false;
        dictionary actuatorModel;
        filterType COrolloffType_;
        scalar errorIntegral_=0;
        scalar errorDifferential_=0.0;
        scalar error_=0;
        scalar controlDelay_;
        bool applyFilter_=false;
        scalar oldError_=0;
        scalar oldErrorIntegral_=0;
        scalar oldOmega_=0.0;
        label  timeIndex_=0;
        scalar oldangle_=0;
        scalar angle=0;
        scalar rawangle_=0;
        scalar rawomega_=0;
        scalar PV_=0; // Process Variable

        int controlTarget_;  // should be label type
        wordList patches_; // did I declare this?
        dictionary forcesDict_;
        scalar setPoint_=-0.35;
        label direction_;
        List<scalar> mycoefs;
        Vector<scalar> myForce; // This is the same as vector class

        //actuator limits
        bool saturate_=false; // 
        scalar anglemax_;
        scalar anglemin_;
        scalar omegamax_;
        scalar omegamin_;
        word controlString_="(Cl)";
        // anti-windup
        scalar satDiff_ = 0;
        bool antiWindup_;
        bool integralReset_;
        //scalar oldSatDiff_ = 0;
        //averaging variables
        scalar meanValue_;
        dictionary PVfilterDict_;
        filterType PVfilterType_;
        scalar window_;
        scalar t0_;
        //scalar Tf_; // 63% response time in first order filter for Control Output
        scalar wc_; //cut-off frequency: widespread concept is where gain=1/sqrt(2) or gain=-3dB
        scalar T_;  //Sampling period for digital filter
        vector COdenom_;//  G(s) = (b0*s^2 + b1*s + b2) / (a0*s^2 + a1*s + a2*s)
        vector COnum_; // if 1st order we will use b1,b2 and a1,a2 and appropriate difference formula
        vector PVdenom_; // MIRA A VER: We could have used the coefficients in G(z) as well
        vector PVnum_;
        vector COinbuffer_; // signal buffers for real-time filtering
        vector COoutbuffer_;
        vector PVinbuffer_;
        vector PVoutbuffer_;
        //RECUERDA HACER EL WARPING
        //functionObject::lforceCoeffs testForces;
    
        
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

        // Test Constructor (for Digital filter test, but customisable for testing other functionalities)
        // Requires a fully defined dynamic mesh case aimed to use this class
        //PIDangularDisplacementPointPatchVectorField(const dictionary&);
        /*PIDangularDisplacementPointPatchVectorField(PIDangularDisplacementPointPatchVectorField& PID):this->clone()
        {
            Info<<"Copy constructor called"<<endl;
        }*/
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
        void rawCO(const scalar dt, const scalar PVsignal);

        //void filterCO(const scalar dt);

        void saturator(const scalar dt);

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            IOdictionary readControl();

            word controlString
            (
                const dictionary&,
                label controlTarget
            );

            void rollavg(scalar currentValue, dictionary& dict, const scalar& dt);

            scalar digitalFilter(scalar inSignal, const vector& num,const vector& denom, vector& inbuf, vector& outbuf, const filterType filterName, scalar& dt);

            tf zplane_tf(vector numerator, vector denominator);
            void selectFilterAndSetup
            (
                vector& num, vector& denom, const filterType filterName, 
                dictionary& filterDict, const dictionary& pdict
            );

        //- Write
        void report(Ostream&) const;
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
