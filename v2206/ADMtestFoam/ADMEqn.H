Info << "ADM convection-diffusion calculation" << endl;



//TODO: turbulent diff and different diff compared to particulate species
volScalarField Deff = turbulence->nu() / reaction->Sc();

{
    // <<< TEST ---------------------------------------------------------------- //

    fvScalarMatrix SubEqn    
    (
        fvm::ddt(Substances)
      + fvm::div(phi, Substances)
      - fvm::laplacian(Deff, Substances)
     ==
        fvOptions(Substances)
    //   + reaction->R(i)
    );

    SubEqn.relax();

    fvOptions.constrain(SubEqn);

    SubEqn.solve(mesh.solver("Yi"));

    fvOptions.correct(Substances);


    // <<< TEST ---------------------------------------------------------------- //

    // forAll(YPtrs, i)
    // {
    //     volScalarField& Yi = YPtrs[i];

    //     fvScalarMatrix YiEqn
    //     (
    //         fvm::ddt(Yi)
    //       + fvm::div(phi, Yi)
    //       - fvm::laplacian(Deff, Yi)
    //      ==
    //         fvOptions(Yi)
    //     //   + reaction->R(i)
    //     );

    //     YiEqn.relax();

    //     fvOptions.constrain(YiEqn);

    //     YiEqn.solve(mesh.solver("Yi"));

    //     fvOptions.correct(Yi);
    // }

    // Gas Phase ============================================================ //
    // TODO: add rho and make adjustment for compressibility

    // forAll(GPtrs, i)
    // {
    //     volScalarField& Gi = GPtrs[i];

    //     fvScalarMatrix GiEqn
    //     (
    //         fvm::ddt(Gi)
    //       + fvm::div(phi, Gi)
    //       - fvm::laplacian(Deff, Gi)
    //      ==
    //         fvOptions(Gi)
    //     //   + reaction->RG(i)
    //     );

    //     GiEqn.relax();

    //     fvOptions.constrain(GiEqn);

    //     GiEqn.solve(mesh.solver("Yi"));

    //     fvOptions.correct(Gi);
    // }
    // Gas Phase ============================================================ //

}