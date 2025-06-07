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
              + fvm::div(testMixture->phiHS(),Yi,"div(phiH,Yi)")
              - fvm::laplacian(testMixture->DalphaS(),Yi)
              ==
                SpY
            );

            // if (phiHScheme == "Gauss linear") YiEqn += fvm::div(phiHS,Yi,"div(phiHS,Yi)");
            // else if (phiHScheme == "Gauss upwind")
            // {
            //     YiEqn += fvm::div(phiHUp,Yi,"div(phiHS,Yi)")
            //     + fvm::div(phiHDown,Yi,"div(phiHS,Yi)");
            // }
            // else
            // {
            //     Info<< "div(phiHS,Yi) should be equal to Gauss linear or Gauss upwind"
            // 	<< endl
            // 	<< abort(FatalError);
            // }

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
}