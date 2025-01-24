{
    if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf))
    );
    MRF.makeRelative(phiHbyA);


    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    if (!correctPhi)
    {
        rAU.clear();
    }
}



// {
//     if (correctPhi)
//     {
//         rAU.ref() = 1.0/UEqn.A();
//     }
//     else
//     {
//         rAU = 1.0/UEqn.A();
//     }

//     surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
//     volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
//     // Info<< ">>> test" << endl;

//     // surfaceScalarField phig(-rAUf*ghf*fvc::snGrad(rhok)*mesh.magSf());
//     surfaceScalarField phig
//     (
//         (
//             mixture.surfaceTensionForce()
//           - ghf*fvc::snGrad(rhok) // <-- rho/rhok ?!!!
//         )*rAUf*mesh.magSf()
//     );


//     surfaceScalarField phiHbyA
//     (
//         "phiHbyA",
//         fvc::flux(HbyA)
//       + MRF.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf))
//     //   + MRF.zeroFilter(rAUf*fvc::ddtCorr(U, phi))
//       + phig
//     );

//     MRF.makeRelative(phiHbyA);

//     // Update the pressure BCs to ensure flux consistency
//     constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

//     while (pimple.correctNonOrthogonal())
//     {
//         fvScalarMatrix p_rghEqn
//         (
//             fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
//         );

//         p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

//         p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

//         if (pimple.finalNonOrthogonalIter())
//         {
//             // Calculate the conservative fluxes
//             phi = phiHbyA - p_rghEqn.flux();

//             // Explicitly relax pressure for momentum corrector
//             p_rgh.relax();

//             // Correct the momentum source with the pressure gradient flux
//             // calculated from the relaxed pressure
//             U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
//             U.correctBoundaryConditions();
//             fvOptions.correct(U);

//             // Correct Uf if the mesh is moving
//             fvc::correctUf(Uf, U, phi);

//             // Make the fluxes relative to the mesh motion
//             fvc::makeRelative(phi, U);
//         }
//     }

//     #include "continuityErrs.H"

//     p = p_rgh + rhok*gh;

//     if (p_rgh.needReference())
//     {
//         p += dimensionedScalar
//         (
//             "p",
//             p.dimensions(),
//             pRefValue - getRefCellValue(p, pRefCell)
//         );
//         p_rgh = p - rhok*gh;
//     }

//     if (!correctPhi)
//     {
//         rAU.clear();
//     }
// }


