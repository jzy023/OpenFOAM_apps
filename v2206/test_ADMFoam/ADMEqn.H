Info<< "ADM convection-diffusion calculation" << endl;

// TODO: implement this with fields tabble in createFields.H
//       to enable divSchemes for Si and Gi [ref: reactingFoam]

// tmp<fv::convectionScheme<scalar>> mvConvection
// (
//     fv::convectionScheme<scalar>::New
//     (
//         mesh,
//         fields,
//         phi,
//         mesh.divScheme("div(phi,Yi_h)")
//     )
// );


//TODO: turbulent diff and different diff compared to particulate species
volScalarField Deff = 
(
    (turbulence->nu() / reaction->Sc())
  + (turbulence->nut() / reaction->Sct())
);


{
    //- Fluid and Particulate Phase
    // const surfaceScalarField alphaRhoPhiLiq(rhoPhi * alphaLiq);
    forAll(YPtrs, i)
    {
        if (i == 7)
        {
            continue;
        }
        
        volScalarField& Yi = YPtrs[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(alphaLiq, rho, Yi)
        //   + fvm::div(rhoPhi, alphaLiq*Yi)
          - fvm::laplacian(alphaLiq*Deff, Yi)
         ==
            fvOptions(Yi)
          + alphaLiq*reaction->R(i)
          + fvc::ddt(/*residualAlpha_*/ 1e-2*rho, Yi)
          - fvm::ddt(/*residualAlpha_*/ 1e-2*rho, Yi)
        );

        YiEqn.relax();

        fvOptions.constrain(YiEqn);

        YiEqn.solve(mesh.solver("Yi"));

        fvOptions.correct(Yi);
    }

    //- Gas Phase
    const volScalarField& alphaGas = 1 - alphaLiq;
    // const surfaceScalarField alphaRhoPhiGas = rhoPhi * alphaGas;
    forAll(GPtrs, i)
    {
        volScalarField& Gi = GPtrs[i];

        fvScalarMatrix GiEqn
        (
            fvm::ddt(alphaGas, Gi)
        //   + fvm::div(rhoPhi, alphaLiq*Gi)
          - fvm::laplacian(alphaGas*Deff, Gi)
         ==
            fvOptions(Gi)
          + alphaGas*reaction->RG(i)
          + fvc::ddt(/*residualAlpha_*/ 1e-2*rho, Gi)
          - fvm::ddt(/*residualAlpha_*/ 1e-2*rho, Gi)
        );

        GiEqn.relax();

        fvOptions.constrain(GiEqn);

        GiEqn.solve(mesh.solver("Yi"));

        fvOptions.correct(Gi);
    }



    //- testing  Gas Phase 
    // forAll(GPtrs_test, i)
    // {
    //     volScalarField& Gi = GPtrs_test[i];

    //     fvScalarMatrix Gi_testEqn
    //     (
    //         fvm::ddt(Gi)
    //       + fvm::div(phi, Gi)
    //       - fvm::laplacian(Deff, Gi)
    //      ==
    //         fvOptions(Gi)
    //       + reaction->RG_test(i)
    //     );

    //     Gi_testEqn.relax();

    //     fvOptions.constrain(Gi_testEqn);

    //     Gi_testEqn.solve(mesh.solver("Yi"));

    //     fvOptions.correct(Gi);
    // }
}