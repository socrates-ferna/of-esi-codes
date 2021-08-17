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

Global
    setMeshDeltaT

Description
    Reset the timestep to maintain a constant maximum Mesh Motion courant Number.
    Reduction of time-step is immediate, but increase is damped to avoid
    unstable oscillations.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    scalar maxMeshDeltaTFact = maxMeshCo/(meshCoNum + SMALL);
    scalar meshdeltaTFact = min(min(maxMeshDeltaTFact, 1.0 + 0.1*maxMeshDeltaTFact), 1.2);


    runTime.setDeltaT
    (
        min
        (
            meshdeltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );

    Info<< "Mesh-corrected deltaT = " <<  runTime.deltaTValue() << endl;
}

// ************************************************************************* //