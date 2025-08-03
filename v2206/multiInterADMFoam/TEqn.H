{
    radiation->correct();
    rhoCp =  rho*mixture.Cp();

    const surfaceScalarField rhoCpPhi(fvc::interpolate(mixture.Cp())*rhoPhi);

    // const volScalarField kappaEff
    // (
    //     "kappaEff",
    //     mixture.kappa() + mixture.Cp()*turbulence->mut()/mixture.Prt()
    // );

    kappaEff = mixture.kappa() + mixture.Cp()*turbulence->mut()/mixture.Prt();

    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T)
      + fvm::div(rhoCpPhi, T, "div(phi,T)")
      - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)
      - fvm::laplacian(kappaEff, T,  "laplacian(kappa,T)")
      ==
        mixture.heatTransfer(T)
      + radiation->ST(T)
      + fvOptions(rhoCp, T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    mixture.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
