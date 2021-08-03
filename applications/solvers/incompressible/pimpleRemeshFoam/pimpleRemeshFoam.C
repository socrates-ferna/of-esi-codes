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
#include "mapFields.H"
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

    IOdictionary myMapDict
    (
        IOobject
        (
            "mymapFieldsDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    //Info<<remeshDict.get<wordRes>("options")<<endl; // crashes bc of -b option
    /*
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
    */
    scalar remeshPeriod(remeshDict.getOrDefault<scalar>("remeshPeriod",0.01));


    
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
        
        if (lastMeshCheck > remeshPeriod)
        {
            mesh.setAspectThreshold(1e10);
            failedChecks = mesh.checkMesh(true);
            if (failedChecks)
            {
                runTime.writeAndEnd();
                
                /*
                Info<<"creating IOobject"<<endl;
                IOobject oMesh
                (
                    "copyMesh",
                    runTime.timeName(), 
                    runTime,
                    Foam::IOobject::NO_READ,
                    Foam::IOobject::AUTO_WRITE
                );

                Info<<"creating auxiliary mesh"<<endl;


                fvMesh auxMesh
                (
                    oMesh,
                    pointField(mesh.points()),
                    faceList(mesh.faces()),
                    cellList(mesh.cells()),
                    true
                );
                
                const polyBoundaryMesh& patches = mesh.boundaryMesh();

                List<polyPatch*> p(patches.size());

                forAll(p,patchi)
                {
                    p[patchi] = patches[patchi].clone(auxMesh.boundaryMesh()).ptr();
                }

                auxMesh.addFvPatches(p);
                
                
                //List<polyPatch*> mypatches = mesh.boundary();
                //auxMesh.addFvPatches(mesh.boundaryMesh(),true);
                auxMesh.write();
                //mesh.resetMotion();
                mesh.write(); // esto me da el problema
                //mesh.destr

                Info<<"Creating mapFields Object with the following dict:"<<endl;
                Info<<myMapDict.subDict("origToCopy")<<endl;
                functionObjects::mapFields myMapFields("oToCopy",runTime,myMapDict.subDict("origToCopy"));

                myMapFields.execute();
                myMapFields.write();
                */
                /*
                mesh.write();
                
                if(Pstream::parRun())
                {
                    if(Pstream::master())
                    {
                        //call remeshing routine
                        system("bash remesh_plus_map");
                        //system("pointwise -b meshUnstructured.glf parameters.dat"); // let's say it writes the mesh to ../aux_case

                                            }
                }
                else
                {
                    system("bash remesh_plus_map");
                    //system("pointwise -b meshUnstructured.glf parameters.dat");
                }
                */
                /*
                //system("pointwise -b meshUnstructured.glf parameters.dat");
                Info<<args.args()<<endl;
                Info<<args.rootPath()<<args.caseName()<<endl;
                //mesh.~dynamicFvMesh();
                //runTime.readIfModified();
                //runTime.readModifiedObjects();
                Info<<"MES SIZE"<<mesh.nCells();
                
                Info<< "RECREATING MESH"<<endl;
                
                autoPtr<dynamicFvMesh> meshPtr(dynamicFvMesh::New(args,runTime));

                dynamicFvMesh& mesh = meshPtr();
                

                */
                //mesh.write();
                //mesh.update();
                
/*
                //#include "createDynamicFvMesh.H"

                functionObjects::mapFields backMapFields("copyMesh",runTime,myMapDict.subDict("copyToNew"));

                backMapFields.execute();
                backMapFields.write();

                auxMesh.~fvMesh();

*/


                // build mesh in aux_case
                // map to that mesh with meshToMesh class
                // "substitute" my current mesh with the mesh built from aux_case
                // map back from aux_case to the current mesh (which is identical to aux_case)



                //reset failedChecks
                //failedChecks = false;
            }
            lastMeshCheck = 0.0;
        }
        
        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
