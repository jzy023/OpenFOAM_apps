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
    forAll(YPtrs, i)
    {
        if (i == 7)
        {
            continue;
        }
        
        volScalarField& Yi = YPtrs[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(Yi)
          + fvm::div(phi, Yi)
          - fvm::laplacian(Deff, Yi)
         ==
            fvOptions(Yi)
          + alphaLiq*reaction->R(i)
        );

        YiEqn.relax();

        fvOptions.constrain(YiEqn);

        YiEqn.solve(mesh.solver("Yi"));

        fvOptions.correct(Yi);
    }

    //- Gas Phase
    forAll(GPtrs, i)
    {
        volScalarField& Gi = GPtrs[i];

        fvScalarMatrix GiEqn
        (
            fvm::ddt(Gi)
          + fvm::div(phi, Gi)
          - fvm::laplacian(Deff, Gi)
         ==
            fvOptions(Gi)
          + alphaGas*reaction->RG(i)
        );

        GiEqn.relax();

        fvOptions.constrain(GiEqn);

        GiEqn.solve(mesh.solver("Yi"));

        fvOptions.correct(Gi);
    }

    // const surfaceScalarField& alphaRhoPhiLiq = phi * fvc::interpolate(alphaLiq * rho);
    // const surfaceScalarField& alphaRhoPhiLiq = rhoPhi * fvc::interpolate(alphaLiq);
    // const surfaceScalarField& alphaRhoPhiLiq = fvc::interpolate(rho)*mixture.phases()["liquid"].alphaPhi();

    //- Fluid and Particulate Phase

    // forAll(YPtrs, i)
    // {
    //     if (i == 7)
    //     {
    //         continue;
    //     }
        
    //     volScalarField& Yi = YPtrs[i];

    //     fvScalarMatrix YiEqn
    //     (
    //         fvm::ddt(alphaLiq, rho, Yi)
    //       + fvm::div(alphaRhoPhiLiq, Yi)
    //       - fvm::laplacian(alphaLiq*rho*Deff, Yi)
    //      ==
    //         alphaLiq*reaction->R(i)
    //     //   + fvOptions(Yi)
    //     //   + fvc::ddt(residualAlpha_*rho, Yi)
    //     //   - fvm::ddt(residualAlpha_*rho, Yi)
    //     );

    //     YiEqn.relax();

    //     fvOptions.constrain(YiEqn);

    //     YiEqn.solve(mesh.solver("Yi"));

    //     fvOptions.correct(Yi);
    // }

    // //- Gas Phase
    
    // const surfaceScalarField alphaRhoPhiGas = rhoPhi * fvc::interpolate(alphaGas);
    
    // forAll(GPtrs, i)
    // {
    //     volScalarField& Gi = GPtrs[i];

    //     fvScalarMatrix GiEqn
    //     (
    //         fvm::ddt(alphaGas, Gi)
    //       + fvm::div(alphaRhoPhiGas, Gi)
    //       - fvm::laplacian(alphaGas*Deff, Gi)
    //      ==
    //         alphaGas*reaction->RG(i)
    //     //   + fvOptions(Gi)
    //     //   + fvc::ddt(residualAlpha_*rho, Gi)
    //     //   - fvm::ddt(residualAlpha_*rho, Gi)
    //     );

    //     GiEqn.relax();

    //     fvOptions.constrain(GiEqn);

    //     GiEqn.solve(mesh.solver("Yi"));

    //     fvOptions.correct(Gi);
    // }

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