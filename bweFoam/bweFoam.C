/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
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
    bweFoam

Description
    Transient solver for Boussinesq-type equations.

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "CourantNo.H"
        #include "setDeltaT.H"


        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;


    //Mass conservation equation
        fvScalarMatrix hEqn(
            fvm::ddt(h)
            + fvm::div(phi,h)
        );

        hEqn.solve();

    // To avoid h<0, we assume a SMALL h when h < hCritical, i.e. H1 = max (H, Hdry)
    //Reference: FUNWAVE-TVD
        H1 = (h-hDry)*pos(h-hDry) + hDry;

	// Bottom friction based on Manning-equation
        fircCoeff = (noDimMagG)*mag(U)*pow(nc,2.0)/pow(H1/dim_m,4.0/3.0)/dim_m;

        wetDry = pos(h-hDry);
    //Momentum conservation equation

        fvVectorMatrix hUEqn
        (
          fvm::ddt(hU)                                                      //Temporal term
	    + fvm::div(phi, hU)                                                 //Spacial term
	    + magg*h*fvc::grad(h + bath)                                        //Flow depth and terrain
	    + fvm::Sp(fircCoeff,hU)                                                 //Bottom friction
	    - fvm::laplacian(Am,hU)                                             //Turbulent stresses
        - h*(1+Beta)*h/2*fvc::grad(fvc::div(h*fvc::ddt(U)))                 //Boussinesq term 1
        - h*Beta*h/2*magg*fvc::grad(fvc::div(h*fvc::grad(h + bath)))        //Boussinesq term 2
        + h*(1+Beta)*h*h/6*fvc::grad(fvc::div(fvc::ddt(U)))                 //Boussinesq term 3
        + h*Beta*h*h/6*magg*fvc::grad(fvc::div(fvc::grad(h + bath)))        //Boussinesq term 4
        );  

        hUEqn.solve();

        hU = hU*pos(H1 - hDry2); 
        U = hU/H1;
        hTotal = h + bath;
       
        phi = (fvc::interpolate(U) & mesh.Sf());

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
