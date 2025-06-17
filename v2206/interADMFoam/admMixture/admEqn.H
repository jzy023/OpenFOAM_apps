{
    Info<< "ADM convection-diffusion calculation" << endl;
    
    // Soluables
    // forAll(Si, i)
    // {
        label i = 6;
        
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
        // Info<< "Species concentration (sum) = "
        // << gSum(Yi.internalField())
        // << "  Min(" << Yi.name() << ") = " << gMin(Yi.internalField())
        // << "  Max(" << Yi.name() << ") = " << gMax(Yi.internalField())
        // << endl;
    // }

    //calculate mass transfer flux
    mixture->limit();
}