    volScalarField rAUc(1.0/UcEqn.A());
    surfaceScalarField rAUcf("Dp", fvc::interpolate(rAUc));
    volVectorField HbyA(constrainHbyA(rAUc*UcEqn.H(), Uc, p));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
           fvc::flux(HbyA)
         // + alphacf*rAUcf*fvc::ddtCorr(Uc, phic)
        )
    );

    if (p.needReference())
    {
        adjustPhi(phiHbyA, Uc, p);
    }

    phiHbyA += phicForces;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, Uc, phiHbyA, rAUcf);

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {   

        fvScalarMatrix pEqn
        (
            fvm::laplacian(alphac*rAUc, p) //**//
         ==
            fvc::ddt(alphac) + fvc::div(alphacf*phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phic = phiHbyA - pEqn.flux()/alphacf;
        }
    }


#include "continuityErrs.H"

//**//
p.relax();

Uc = HbyA
   + rAUc*fvc::reconstruct(phicForces/rAUcf)
   - rAUc*fvc::grad(p)/alphac;

Uc.correctBoundaryConditions();
