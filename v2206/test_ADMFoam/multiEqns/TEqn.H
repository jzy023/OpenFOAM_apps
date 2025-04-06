{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    // TODO: Cp or both phases?
    volScalarField alphaEff("alphaEff", mixture.nu()/Pr + alphat);
    volScalarField kappaEff("kappaEff", rho*alphaEff);

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T, "div(phi,T)")
      - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), T)
      - fvm::laplacian(kappaEff, T)
      ==
        fvOptions(rho, T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    mixture.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;

    // rhok = 1.0 - beta*(T - TRef);
}

// TODO:
// consult interCondensatingEvaporatingFoam/TEqn.H