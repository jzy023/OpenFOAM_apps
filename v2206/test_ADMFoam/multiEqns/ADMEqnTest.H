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
    // volScalarField& Yi = YPtrs[6];
    // for (int i = 0; i < 3; i++)
    // {
    //     fvScalarMatrix YiEqn
    //     (
    //         fvm::ddt(Yi)
    //       + fvm::div(phi, Yi, "div(phi,test)")
    //       - fvm::laplacian(Deff, Yi)
    //     ==
    //         alphaLiq*fvOptions(Yi)
    //     );

    //     YiEqn.relax();
    //     fvOptions.constrain(YiEqn);
    //     YiEqn.solve();
    //     fvOptions.correct(Yi);

    //     tPhiYi = YiEqn.flux();
    // }

    // MULES::explicitSolve
    // (
    //     geometricOneField(),
    //     Yi,
    //     phi,
    //     tPhiYi.ref(),
    //     oneField(),
    //     zeroField()
    // );

    // --------------------------------------------------------------------------------------------------

        volScalarField& Yi = YPtrs[6];

        // Deff *= pos(alphaLiq - 0.99);
        // tmp<surfaceScalarField> tPhiYi;

        for (int i = 0; i < 2; i++)
        {
            fvScalarMatrix YiEqn
            (
                fvm::ddt(Yi)
              - fvm::laplacian(alphaLiq*Deff, Yi)
              + fvm::div(mixture.phases()["liquid"].alphaPhi(), Yi, "div(phi,test)")
            //   + fvm::div(phi, Yi, "div(phi,test)")
            ==
                alphaLiq*fvOptions(Yi)
            );

            YiEqn.relax();
            fvOptions.constrain(YiEqn);
            YiEqn.solve();
            fvOptions.correct(Yi);

            // tPhiYi = YiEqn.flux();
        }

        // MULES::explicitSolve
        // (
        //     geometricOneField(),
        //     Yi,
        //     phi,
        //     tPhiYi.ref(),
        //     oneField(),
        //     zeroField()
        // );

    // --------------------------------------------------------------------------------------------------

        //- Fluid and Particulate Phase

        // const volScalarField residualAlpha_
        // (
        //     IOobject
        //     (
        //         "residualAlpha",
        //         mesh.time().timeName(),
        //         mesh,
        //         IOobject::NO_READ, // READ_IF_PRESENT,
        //         IOobject::NO_WRITE
        //     ),
        //     mesh,
        //     dimensionedScalar
        //     (
        //         dimless,
        //         1e-3
        //     )
        // );
            
        // volScalarField& Yi = YPtrs[6];

        // // Deff *= pos(alphaLiq - 0.99);

        // for (size_t i = 0; i < 3; i++)
        // {
        //     fvScalarMatrix YiEqn
        //     (
        //         fvm::ddt(alphaLiq, Yi)
        //       - fvm::laplacian(alphaLiq*Deff, Yi)
        //       + fvm::div(mixture.phases()["liquid"].alphaPhi(), Yi, "div(phi,test)")
        //     //   + fvm::div(phi, Yi, "div(phi,test)")
        //     ==
        //     //     alphaLiq*fvOptions(Yi)
        //        fvc::ddt(residualAlpha_, Yi)
        //      - fvm::ddt(residualAlpha_, Yi)
        //     );

        //     YiEqn.relax();
        //     fvOptions.constrain(YiEqn);
        //     YiEqn.solve();
        //     fvOptions.correct(Yi);
        // }


    // --------------------------------------------------------------------------------------------------

        // - Fluid and Particulate Phase
        
        // volScalarField& Yi = YPtrs[6];
        // volScalarField YiAlpha = Yi;
        // YiAlpha.field() = max(alphaLiq, scalar(1e-3))().field() * Yi.field();

        // Deff *= pos(alphaLiq - 0.99);

        // for (size_t i = 0; i < 3; i++)
        // {
        //     fvScalarMatrix YiEqn
        //     (
        //         fvm::ddt(YiAlpha)
        //     //   + fvm::div(phi, YiAlpha, "div(phi,test)")
        //       + fvm::div(mixture.phases()["liquid"].alphaPhi(), Yi, "div(phi,test)")
        //       - fvm::laplacian(Deff, YiAlpha)
        //     ==
        //         alphaLiq*fvOptions(YiAlpha)
        //     );

        //     YiEqn.relax();
        //     fvOptions.constrain(YiEqn);
        //     YiEqn.solve();
        //     fvOptions.correct(YiAlpha);
        //     Yi = YiAlpha / max(alphaLiq, scalar(1e-3));
        // }
    

    // --------------------------------------------------------------------------------------------------

    // const volScalarField residualAlpha_
    // (
    //     IOobject
    //     (
    //         "residualAlpha",
    //         mesh.time().timeName(),
    //         mesh,
    //         IOobject::NO_READ, // READ_IF_PRESENT,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar
    //     (
    //         dimless,
    //         1e-3
    //     )
    // );

    // const volScalarField rhoLiq
    // (
    //     IOobject
    //     (
    //         "rhoLiq",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     mixture.phases()["liquid"].rho()
    // );

    // const surfaceScalarField& alphaRhoPhiLiq = fvc::interpolate(rhoLiq)*mixture.phases()["liquid"].alphaPhi();

    // //- Fluid and Particulate Phase
        
    //     volScalarField& Yi = YPtrs[6];

    //     fvScalarMatrix YiEqn
    //     (
    //         fvm::ddt(alphaLiq, rhoLiq, Yi)
    //       + fvm::div(alphaRhoPhiLiq, Yi, "div(alphaRhoPhi.*,Yi)")
    //       - fvm::laplacian
    //         (
    //             fvc::interpolate(alphaLiq)*
    //             fvc::interpolate(rhoLiq*Deff),
    //             Yi
    //         )
    //      ==
    //         fvc::ddt(residualAlpha_*rhoLiq, Yi)
    //       - fvm::ddt(residualAlpha_*rhoLiq, Yi)
    //     );

    //     YiEqn.relax();
    //     fvOptions.constrain(YiEqn);
    //     YiEqn.solve();
    //     fvOptions.correct(Yi);

}