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
            reaction->R(i)
          + fvOptions(Yi)
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
            reaction->RG(i)
          + fvOptions(Gi)
        );

        GiEqn.relax();

        fvOptions.constrain(GiEqn);

        GiEqn.solve(mesh.solver("Yi"));

        fvOptions.correct(Gi);
    }
}