/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Author: Socrates Fernandez
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

Application
    pimpleRemeshFoam.C

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    on a moving mesh.

    \heading Solver details
    The solver uses the PIMPLE (merged PISO-SIMPLE) algorithm to solve the
    continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}} + \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

Note
   The motion frequency of this solver can be influenced by the presence
   of "updateControl" and "updateInterval" in the dynamicMeshDict.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
//#include "IOdictionary.H"
//#include "argList.H"

//#include "checkTools.H"
//#include "checkTopology.H"
//#include "checkGeometry.H"
//#include "checkMeshQuality.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " of Newtonian fluids on a moving mesh."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
/*
    if (args.found("parallel"))
    {
        label n = 0;
        //label n = Pstream::nProcs(worldComm);
        Info << "Running on parallel on " << n << "processes"<< endl;
    }

    if (Pstream::parRun())
    {
        Info << "esta forma tb funciona"<< endl;
    }
    if (Pstream::master())
    {
        Info << "I am the master process"<< endl;
    }
*/
    turbulence->validate();

    //label nFailedChecks = 0;

    bool failedChecks = false;

    scalar maxMeshCo = 1.0;

    scalar meshCoNum = 1.0;

    scalar lastMeshCheck = 0.0;

    // --- Read mesh quality criteria

    Info<< "\nReading Mesh Quality Dict\n"<< endl;
    autoPtr<IOdictionary> qualDict;
    //autoPtr<surfaceWriter> surfWriter;
    qualDict.reset
    (
        new IOdictionary
        (
            IOobject
            (
                "meshQualityDict",
                mesh.time().system(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    IOdictionary remeshDict
    (
        IOobject
        (
            "remeshDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    //Info<<remeshDict.get<wordRes>("options")<<endl; // crashes bc of -b option
    string remeshCommand;
    string mesher(remeshDict.getOrDefault<word>("mesher","snappyHexMesh"));
    string language(remeshDict.getOrDefault<word>("language","bash"));
    string logFile(remeshDict.getOrDefault<word>("logFile","log.pimpleRemeshFoam"));
    string script(remeshDict.getOrDefault<word>("script","remesh.sh"));
    string op1(remeshDict.getOrDefault<word>("option1",""));
    op1 = "-" + op1;
    string op2(remeshDict.getOrDefault<word>("option2",""));
    string op3(remeshDict.getOrDefault<word>("option3",""));
    remeshCommand = language+ " " + script+ " " + logFile+ " " + mesher+ " " + op1+ " " + op2+ " " + op3;
    Info<<"Remesh command is: "<<remeshCommand<<endl;
    Info<< "Test call"<<endl;
    system(remeshCommand);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"   //MIRA ESTO POR SI TE ES ÚTIL
        #include "CourantNo.H"
        #include "setDeltaT.H"
        #include "setMeshDeltaT.H"
        ++runTime;
        lastMeshCheck = lastMeshCheck + runTime.deltaTValue();
        Info<< "Time = " << runTime.timeName() << nl << endl; 
        //Info<< "Current time index is: "<< runTime.timeIndex()<< nl << endl;
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }
        //nFailedChecks = checkMeshQuality(mesh, qualDict(), surfWriter);
        //mesh check based on 2 conditions:
        //1. Every 0.1s
        if (lastMeshCheck > 0.1)
        {
            failedChecks = mesh.checkMesh(true);
            if (failedChecks)
            {
                if(Pstream::parRun())
                {
                    if(Pstream::master())
                    {
                        //call remeshing routine
                        system(remeshCommand);
                
                        //call mapFields
                        system("mapFields -case ../aux_case -consistent -sourceTime latestTime .");
                        //missing angle position for correct remeshing and mapFieldsDict for cuttingPatches
                    }
                }
                else
                {
                    system(remeshCommand);
                    system("mapFields -case ../aux_case -consistent -sourceTime latestTime .");
                }
                //reset failedChecks
                failedChecks = false;
            }
            lastMeshCheck = 0.0;
        }
        
        //determine whether remesh is necesary
        //build a mapPolyMesh object
        //remesh
        //mapFields
        /*Info<< "\n Failed " << failedChecks << " mesh checks. \n"<< endl;
        if (failedChecks > 0)
        {
            Info<< "\nCalling remesh routine.\n"<< endl;
            //system("pointwise -b meshUnstructured.glf parameters.dat");
        }
        */
        //mesh.readIfModified();
        //runTime.lookupObject()
        
        //system("mapFields -case ../aux_case -consistent -sourceTime latestTime .");
        
        //system("mapFields -consistent -sourceTime latestTime ../aux_case");
        //runTime.readIfModified();
        //const IOdictionary& dynamicMeshDict = mesh.lookupObject<IOdictionary>("dynamicMeshDict");
        //Info<< "NAMES "<< mesh.names()<< endl;
        //Info<< "path: "<< dynamicMeshDict.filePath()<< endl;
        //Info<< "dynamicMeshDict modified: "<< dynamicMeshDict.modified()<< endl;
        //static_cast<Foam::IOdictionary>(pdisp).write(Info,true);
        //word myname = "forceCoeffs";
        //pdisp.writeData(Info);
        //functionObjectList runTime.functionObjects(myList);
        //Info<< myList.get << endl;
        //bool test = runTime.foundObject<functionObject>(myname);
        //Info<< test <<endl;
        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
