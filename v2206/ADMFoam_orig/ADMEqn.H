Info<< "ADM convection-diffusion calculation" << endl;



//TODO: turbulent diff and different diff compared to particulate species
// volScalarField DpartEff = Deff / 10.0;
volScalarField Deff = 
(
    (turbulence->nu() / reaction->Sc())
  + (turbulence->nut() / reaction->Sct())
);

{
    // <<< TEST ---------------------------------------------------------------- //
    // volScalarField rho
    // (
    //     IOobject
    //     (
    //         "rho",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::READ_IF_PRESENT,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar
    //     (
    //         "rho", 
    //         dimensionSet(1,-3,0,0,0,0,0), 
    //         1000
    //     )
    // );
    // volScalarField DMuEff = Deff*rho;
    // <<< TEST ---------------------------------------------------------------- //

    //- TEST
    // volScalarField test(fvc::div(phi, 1000.0*YPtrs[0]));
    // Info<< test << endl;

    // Qdot = reaction->Qdot();

    // volScalarField& Yi = YPtrs[7];
    // Info<< Yi.dimensions() << endl;
    // Info<< fvm::ddt(Yi)->dimensions() << endl;
    // Info<< fvm::div(phi, Yi)->dimensions() << endl;
    // Info<< reaction->R(7)->dimensions() << endl;

    // reaction->clear();
    // reaction->correct(phi, T);
    // PtrList<volScalarField>& YPtrs = reaction->Y();
    // PtrList<volScalarField>& GPtrs = reaction->G();
    // PtrList<volScalarField>& IOPtrs = reaction->IO();

    forAll(YPtrs, i)
    {
        if (i == 7)
        {
            continue;
        }
        
        // if (i != 7) 
        // {
            volScalarField& Yi = YPtrs[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(Yi)
              + fvm::div(phi, Yi)
              - fvm::laplacian(Deff, Yi)
             ==
                fvOptions(Yi)
              + reaction->R(i)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);
        // }
        // else
        // {
        //     volScalarField& Yi = YPtrs[i];

        //     fvScalarMatrix YiEqn
        //     (
        //         fvm::ddt(Yi)
        //       + fvm::div(phi, Yi)
        //       - fvm::laplacian(DpartEff, Yi)
        //      ==
        //         fvOptions(Yi)
        //       + reaction->R(i)
        //     );

        //     YiEqn.relax();

        //     fvOptions.constrain(YiEqn);

        //     YiEqn.solve(mesh.solver("Yi"));

        //     fvOptions.correct(Yi);
        // }
    }

    // Gas Phase ============================================================ //
    // TODO: add rho and make adjustment for compressibility

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
          + reaction->RG(i)
        );

        GiEqn.relax();

        fvOptions.constrain(GiEqn);

        GiEqn.solve(mesh.solver("Yi"));

        fvOptions.correct(Gi);
    }

    // Ions ================================================================= //
    // TODO: add rho and make adjustment for compressibility

    // forAll(IOPtrs, i)
    // {
    //     volScalarField& IOi = IOPtrs[i];

    //     fvScalarMatrix IOiEqn
    //     (
    //         fvm::ddt(IOi)
    //      ==
    //         reaction->RIO(i)
    //     );

    //     IOiEqn.relax();

    //     fvOptions.constrain(IOiEqn);

    //     IOiEqn.solve(mesh.solver("Yi"));

    //     fvOptions.correct(IOi);
    // }

    // Gas Phase ============================================================ //

}