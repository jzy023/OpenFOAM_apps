{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);
    volScalarField kappaEff("kappaEff", rho*alphaEff);

    // fvScalarMatrix TEqn
    // (
    //     fvm::ddt(T)
    //   + fvm::div(phi, T)
    //   - fvm::laplacian(alphaEff, T)
    //  ==
    //     fvOptions(T)
    //   + radiation->ST(rhoCpRef, T)
    // );

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T, "div(phi,T)")
      - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), T)
      - fvm::laplacian(kappaEff, T)
      ==
        fvOptions(rho, T)
    //   + fluid.heatTransfer(T)
    //   + radiation->ST(T)
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