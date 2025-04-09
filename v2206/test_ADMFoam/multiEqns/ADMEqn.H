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
    
    // TODO: add dict entry?
    int nCorr = 1;

    forAll(YPtrs, i)
    {
        if (YPtrs[i].name() == "Sh2")
        {
            continue;
        }

        volScalarField& Yi = YPtrs[i];

        for (int i = 0; i < nCorr; i++)
        {
            fvScalarMatrix YiEqn
            (
                fvm::ddt(Yi)
              - fvm::laplacian(alphaLiq*Deff, Yi)
              + fvm::div(mixture.phases()["liquid"].alphaPhi(), Yi, "div(phi,Yi)")
            ==
                alphaLiq*fvOptions(Yi)
            );

            YiEqn.relax();
            fvOptions.constrain(YiEqn);
            YiEqn.solve();
            fvOptions.correct(Yi);
        }
    }

    //- Gas Phase
    
    forAll(GPtrs, i)
    {
        volScalarField& Gi = GPtrs[i];

        for (int i = 0; i < nCorr; i++)
        {
            fvScalarMatrix GiEqn
            (
                fvm::ddt(Gi)
              - fvm::laplacian(alphaGas*Deff, Gi)
              + fvm::div(mixture.phases()["gas"].alphaPhi(), Gi, "div(phi,Gi)")
            ==
                alphaGas*fvOptions(Gi)
            );

            GiEqn.relax();
            fvOptions.constrain(GiEqn);
            GiEqn.solve();
            fvOptions.correct(Gi);
        }
    }
}