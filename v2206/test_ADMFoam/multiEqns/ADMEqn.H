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

    const volScalarField residualAlpha_
    (
        IOobject
        (
            "residualAlpha",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ, // READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            dimless,
            1e-3
        )
    );

    // forAll(YPtrs, i)
    // {
    //     if (i == 7)
    //     {
    //         continue;
    //     }
        
    //     volScalarField& Yi = YPtrs[i];

    //     fvScalarMatrix YiEqn
    //     (
    //         fvm::ddt(alphaLiq, Yi)
    //       + fvm::div(mixture.phases()["liquid"].alphaPhi(), Yi)
    //       - fvm::laplacian(alphaLiq*Deff, Yi)
    //      ==
    //         // alphaLiq*reaction->R(i)
    //         alphaLiq*fvOptions(Yi)
    //       + fvc::ddt(residualAlpha_, Yi)
    //       - fvm::ddt(residualAlpha_, Yi)
    //     );

    //     YiEqn.relax();

    //     fvOptions.constrain(YiEqn);

    //     YiEqn.solve(mesh.solver("Yi"));

    //     fvOptions.correct(Yi);
    // }


    // //- Gas Phase
    
    // forAll(GPtrs, i)
    // {
    //     volScalarField& Gi = GPtrs[i];

    //     fvScalarMatrix GiEqn
    //     (
    //         fvm::ddt(alphaGas, Gi)
    //       + fvm::div(mixture.phases()["gas"].alphaPhi(), Gi)
    //       - fvm::laplacian(alphaGas*Deff, Gi)
    //      ==
    //         // alphaGas*reaction->RG(i)
    //         alphaGas*fvOptions(Gi)
    //       + fvc::ddt(residualAlpha_, Gi)
    //       - fvm::ddt(residualAlpha_, Gi)
    //     );

    //     GiEqn.relax();

    //     fvOptions.constrain(GiEqn);

    //     GiEqn.solve(mesh.solver("Gi"));

    //     fvOptions.correct(Gi);
    // }


    // //- Fluid and Particulate Phase
    // forAll(YPtrs, i)
    // {
    //     if (i == 7)
    //     {
    //         continue;
    //     }
        
    //     volScalarField& Yi = YPtrs[i];
    //     volScalarField YiAlpha = Yi;
    //     YiAlpha.field() = max(alphaLiq, scalar(1e-3))().field() * Yi.field();

    //     fvScalarMatrix YiEqn
    //     (
    //         fvm::ddt(YiAlpha)
    //       + fvm::div(phi, YiAlpha)
    //       - fvm::laplacian(Deff, YiAlpha)
    //      ==
    //     //   + alphaLiq*reaction->R(i)
    //         fvOptions(YiAlpha)
    //     );

    //     YiEqn.relax();

    //     fvOptions.constrain(YiEqn);

    //     YiEqn.solve(mesh.solver("Yi"));

    //     fvOptions.correct(YiAlpha);

    //     Yi = YiAlpha / max(alphaLiq, scalar(1e-3));
    // }

    // //- Gas Phase
    // forAll(GPtrs, i)
    // {
    //     volScalarField& Gi = GPtrs[i];
    //     volScalarField GiAlpha = Gi;
    //     GiAlpha.field() = max(alphaGas, scalar(1e-3))().field() * Gi.field();

    //     fvScalarMatrix GiEqn
    //     (
    //         fvm::ddt(GiAlpha)
    //       + fvm::div(phi, GiAlpha)
    //       - fvm::laplacian(Deff, GiAlpha)
    //      ==
    //     //   + alphaGas*reaction->RG(i)
    //         fvOptions(GiAlpha)
    //     );

    //     GiEqn.relax();

    //     fvOptions.constrain(GiEqn);

    //     GiEqn.solve(mesh.solver("Yi"));

    //     fvOptions.correct(GiAlpha);

    //     Gi = GiAlpha / max(alphaGas, scalar(1e-3));
    // }

    // const surfaceScalarField& alphaRhoPhiLiq = phi * fvc::interpolate(alphaLiq * rho);
    // const surfaceScalarField& alphaRhoPhiLiq = rhoPhi * fvc::interpolate(alphaLiq);
    // const dimensionedScalar& rhoLiq = mixture.phases()["liquid"].rho();
    // const dimensionedScalar& rhoGas = mixture.phases()["gas"].rho();

    const volScalarField rhoLiq
    (
        IOobject
        (
            "rhoLiq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        mixture.phases()["liquid"].rho()
    );

    const volScalarField rhoGas
    (
        IOobject
        (
            "rhoGas",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        mixture.phases()["gas"].rho()
    );

    const surfaceScalarField& alphaRhoPhiLiq = fvc::interpolate(rhoLiq)*mixture.phases()["liquid"].alphaPhi();
    const surfaceScalarField& alphaRhoPhiGas = fvc::interpolate(rhoGas)*mixture.phases()["gas"].alphaPhi();

    // const surfaceScalarField& alphaRhoPhiLiq = fvc::interpolate(rhoLiq)*mixture.phases()["liquid"].alphaPhi();
    // const surfaceScalarField& alphaRhoPhiGas = fvc::interpolate(rhoGas)*mixture.phases()["gas"].alphaPhi();


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
            fvm::ddt(alphaLiq, rhoLiq, Yi)
          + fvm::div(alphaRhoPhiLiq, Yi, "div(alphaRhoPhi.*,Yi)")
          - fvm::laplacian
            (
                fvc::interpolate(alphaLiq)*
                fvc::interpolate(rhoLiq*Deff),
                Yi
            )
         ==
            // alphaLiq*rhoLiq*reaction->R(i)
        //   + fvOptions(Yi)
            fvc::ddt(residualAlpha_*rhoLiq, Yi)
          - fvm::ddt(residualAlpha_*rhoLiq, Yi)
        );

        YiEqn.relax();

        fvOptions.constrain(YiEqn);

        YiEqn.solve(mesh.solver("Yi"));

        fvOptions.correct(Yi);
    }

    // //- Gas Phase
    
    // forAll(GPtrs, i)
    // {
    //     volScalarField& Gi = GPtrs[i];

    //     fvScalarMatrix GiEqn
    //     (
    //         fvm::ddt(alphaGas, rhoGas, Gi)
    //       + fvm::div(alphaRhoPhiGas, Gi)
    //       - fvm::laplacian
    //         (
    //             fvc::interpolate(alphaGas)*
    //             fvc::interpolate(rhoGas*Deff),
    //             Gi
    //         )
    //      ==
    //     //     alphaGas*rhoGas*reaction->RG(i)
    //     //   + fvOptions(Gi)
    //         fvc::ddt(residualAlpha_*rhoGas, Gi)
    //       - fvm::ddt(residualAlpha_*rhoGas, Gi)
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