{
    Info<< "ADM convection-diffusion calculation" << endl;
    
    // Soluables
    forAll(Si, i)
    {
        volScalarField& Yi = Si[i];

        volScalarField SpY = (Yi - Yi.oldTime()) / runTime.deltaT();

        while (pimple.correctNonOrthogonal())
        {
            //compute Y equation
            fvScalarMatrix YiEqn
            (
                fvm::ddt(Yi)
              + fvm::div(mixture->phiHS(),Yi,"div(phiH,Yi)")
              - fvm::laplacian(mixture->DalphaS(),Yi)
             ==
                SpY
            );

            //solve equations
            if (pimple.finalNonOrthogonalIter())
            {
                YiEqn.solve(mesh.solver("YiFinal"));
            }
            else
            {
                YiEqn.solve(mesh.solver("Yi"));
            }
        }

        // Info<< "Species concentration = "
        // << Yi.weightedAverage(mesh.V()).value()
        Info<< "Species concentration (sum) = "
        << gSum(Yi.internalField())
        << "  Min(" << Yi.name() << ") = " << gMin(Yi.internalField())
        << "  Max(" << Yi.name() << ") = " << gMax(Yi.internalField())
        << endl;
    }

    // testing ------------------------------------------------------------
    // TODO: maybe make a function mixture->calcReactionGas(runTime.deltaT());
    //                                                      dimensionedScalar
    PtrList<volScalarField::Internal>& dG = mixture->reaction()->dG_test();
    forAll(dG, i)
    {
        mixture->reaction()->G_test()[i] += 
        (
            dG[i].weightedAverage(mesh.V()).value() * runTime.deltaT()
        );
    }

    // //calculate mass transfer flux
    // mixture->limit();

    // // Gaseous
    // // forAll(Gi, i)
    // // {
    //     label i = 0;

    //     volScalarField& Yi = Gi[i];

    //     volScalarField SpY = (Yi - Yi.oldTime()) / runTime.deltaT();

    //     while (pimple.correctNonOrthogonal())
    //     {
    //         //compute Y equation
    //         fvScalarMatrix YiEqn
    //         (
    //             fvm::ddt(Yi)
    //           + fvm::div(mixture->phiHS(),Yi,"div(phiH,Yi)")
    //           - fvm::laplacian(mixture->DalphaS(),Yi)
    //          ==
    //             SpY
    //         );

    //         //solve equations
    //         if (pimple.finalNonOrthogonalIter())
    //         {
    //             YiEqn.solve(mesh.solver("YiFinal"));
    //         }
    //         else
    //         {
    //             YiEqn.solve(mesh.solver("Yi"));
    //         }
    //     }

    //     // Info<< "Species concentratio n = "
    //     // << Yi.weightedAverage(mesh.V()).value()
    //     // Info<< "Species concentration (sum) = "
    //     // << gSum(Yi.internalField())
    //     // << "  Min(" << Yi.name() << ") = " << gMin(Yi.internalField())
    //     // << "  Max(" << Yi.name() << ") = " << gMax(Yi.internalField())
    //     // << endl;
    // // }
}