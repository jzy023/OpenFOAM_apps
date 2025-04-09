/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |  
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C)  Jeremy Z. Yan
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
    ADMFoam.C

Group
    grpIncompressibleSolvers

Description
    Transient solver for ADM no1 coupled with incompressible, turbulent flow of 
    Newtonian fluids on a moving mesh.

    \heading Solver details
    The solver uses the PIMPLE (merged PISO-SIMPLE) algorithm to solve the
    continuity equation

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

Note

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvcSmooth.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "radiationModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"

#include "localEulerDdtScheme.H"
#include "ADMno1.H"

// testing
// > multiphaseInterFoam
#include "multiphaseADMixture.H"
#include "CMULES.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for ADM no1 coupled with incompressible," 
        " turbulent flow of Newtonian fluids on a moving mesh."
    );

    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    
    // #include "createControl.H"
    // #include "createTimeControls.H"

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // testing
    // > multiphaseInterFoam
    //- some references to multiphase simulation
    #include "initCorrectPhi.H"
    const surfaceScalarField& rhoPhi(mixture.rhoPhi());
    const volScalarField& alphaLiq = mixture.phases()["liquid"];
    const volScalarField& alphaGas = mixture.phases()["gas"];



    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // ADM1 reaction source terms
    PtrList<volScalarField>& YPtrs = reaction->Y(); // <-- TODO: const?
    PtrList<volScalarField>& GPtrs = reaction->G(); // <-- TODO: const?

    // Loop
    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"


        // testing
        // > multiphaseInterFoam
        #include "alphaCourantNo.H"
        

        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // ADM1 reaction source terms
        reaction->clear();
        reaction->correct
        (   // <-- phiAlpha? or add "alpha"
            phi, 
            alphaLiq,
            alphaGas,
            T,
            p
        ); 
        // reaction->correct(phi, T); // <-- phiAlpha? or add "alpha"


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
            
            // testing
            mixture.solve(reaction->vDotList_test);
            rho = mixture.rho();
            
            #include "multiEqns/UEqn.H"

            #include "multiEqns/TEqn.H"
            
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "multiEqns/pEqn.H"
            }

            // --- ADM calculation
            #include "multiEqns/ADMEqn.H"
            // #include "multiEqns/ADMEqnTest.H"

            if (pimple.turbCorr())
            {
                // laminarTransport.correct();
                // testing
                turbulence->correct();
            }

        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
