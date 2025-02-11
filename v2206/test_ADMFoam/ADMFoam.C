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

// > icoReactingMultiphaseInterFoam
// #include "multiphaseSystem.H"
// #include "turbulentFluidThermoModel.H"


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
    #include "initCorrectPhi.H"
    const surfaceScalarField& rhoPhi(mixture.rhoPhi());

    // > icoReactingMultiphaseInterFoam
    // #include "createFieldRefs.H"
    // #include "createFvOptions.H"



    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"


        // testing
        // > multiphaseInterFoam
        #include "alphaCourantNo.H"

        // > icoReactingMultiphaseInterFoam
        // #include "icoReactingMixture/alphaCourantNo.H"


        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

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

            // #include "TEqn.H"     
            
            // testing
            // > multiphaseInterFoam
            mixture.solve(reaction->vDotList_test);
            rho = mixture.rho();

            // #include "alphaControls.H"      // <-- !!!
            // #include "alphaEqnSubCycle.H"   // <-- !!!

            // > icoReactingMultiphaseInterFoam
            


            // #include "UEqn.H"
            #include "multiphaseADMixture/UEqn.H"
            // #include "icoReactingMixture/UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                // #include "pEqn.H"
                #include "multiphaseADMixture/pEqn.H"
                // #include "icoReactingMixture/pEqn.H"
            }

            if (pimple.turbCorr())
            {
                // laminarTransport.correct();
                // testing
                turbulence->correct();
            }

        }

        // ADM1 reaction source terms
        reaction->clear();
        reaction->correct(phi, T, p);
        // reaction->correct(phi, T);

        PtrList<volScalarField>& YPtrs = reaction->Y();
        PtrList<volScalarField>& GPtrs = reaction->G();
        // PtrList<volScalarField>& GPtrs_test = reaction->G_test();
        
        // --- ADM calculation
        #include "ADMEqn.H"

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
