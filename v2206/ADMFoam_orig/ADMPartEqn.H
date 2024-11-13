Info<< "ADM convection-diffusion calculation for particulates" << endl;

//TODO: turbulent diff and different diff compared to soluable species
// volScalarField Deff = turbulence->nu() / reaction->Sc();
// volScalarField DPartEff = Deff / 10.0;

{
    for (label i = 12; i < 24; i++)
    {
        volScalarField& Yi = YPtrs[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(Yi)
          + fvm::div(phi, Yi)
        //   - fvm::laplacian(DPartEff, Yi)
         ==
            fvOptions(Yi)
        );

        YiEqn.relax();

        fvOptions.constrain(YiEqn);

        YiEqn.solve(mesh.solver("Yi"));

        fvOptions.correct(Yi);
    }       
}