Info << "ADM convection-diffusion calculation" << endl;

PtrList<volScalarField>& YPtrs = reaction->Y();
PtrList<volScalarField>& GPtrs = reaction->G();

//TODO: turbulent diff and different diff compared to particulate species
volScalarField Deff = turbulence->nu() / reaction->Sc();
volScalarField DpartEff = Deff / 10.0;

{
    // Qdot = reaction->Qdot();

    // volScalarField& Yi = YPtrs[7];
    // Info << Yi.dimensions() << endl;
    // Info << fvm::ddt(Yi)->dimensions() << endl;
    // Info << fvm::div(phi, Yi)->dimensions() << endl;
    // Info << reaction->R(7)->dimensions() << endl;

    // reaction->clear();
    // reaction->correct(phi, T);

    forAll(YPtrs, i)
    {
            volScalarField& Yi = YPtrs[i];

            fvScalarMatrix YiEqn
            (
                fvm::div(phi, Yi)
              - fvm::laplacian(Deff, Yi)
             ==
                fvOptions(Yi)
            //   + reaction->R(i)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);
    }

    // Gas Phase ============================================================ //
    // TODO: add rho and make adjustment for compressibility

    forAll(GPtrs, i)
    {
        volScalarField& Gi = GPtrs[i];

        fvScalarMatrix GiEqn
        (
            fvm::div(phi, Gi)
          - fvm::laplacian(Deff, Gi)
         ==
            fvOptions(Gi)
        //   + reaction->RG(i)
        );

        GiEqn.relax();

        fvOptions.constrain(GiEqn);

        GiEqn.solve(mesh.solver("Yi"));

        fvOptions.correct(Gi);
    }

}