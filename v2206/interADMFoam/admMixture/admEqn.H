{
    Info<< "ADM convection-diffusion calculation" << endl;

    volScalarField limitedAlpha1
    (
        min(max(alpha1, scalar(0)), scalar(1))
    );
    
    // Soluables
    PtrList<volScalarField>& SiAlpha = mixture->SiAlpha(); // <-- TODO: const?

    forAll(SiAlpha, i)
    {
        // TODO: simplify?
        volScalarField& YiAlpha = SiAlpha[i];

        if (i == 7)
        {   // TODO: maybe clean up? 
            YiAlpha = limitedAlpha1*mixture->reaction()->Y()[i];
            YiAlpha.correctBoundaryConditions();
            continue;
        }

        volScalarField SpY = (YiAlpha - YiAlpha.oldTime()) / runTime.deltaT();

        while (pimple.correctNonOrthogonal())
        {
            //compute Y equation
            fvScalarMatrix YiEqn
            (
                fvm::ddt(YiAlpha)
              + fvm::div(mixture->phiHS(),YiAlpha,"div(phiH,Yi)")
              - fvm::laplacian(mixture->DalphaS(),YiAlpha)
             ==
                limitedAlpha1*mixture->dSi()[i]
              + SpY
            );

            //solve equations
            if (pimple.finalNonOrthogonalIter())
            {
                YiEqn.solve(mesh.solver("YiFinal"));
            }
            else
            {
                YiEqn.solve(mesh.solver("YiAlpha"));
            }
        }
    }

    // Gaseous
    mixture->updateReactionGas(runTime.deltaT());

    // DEBUG
    // forAll(SiAlpha, i)
    // {
    //     // TODO: simplify?
    //     volScalarField& YiAlpha = SiAlpha[i];

    //     Info<< YiAlpha.name() << " concentration (sum) = "
    //         << gSum(YiAlpha.internalField())
    //         << "  Min = " << gMin(YiAlpha.internalField())
    //         << "  Max = " << gMax(YiAlpha.internalField())
    //         << endl;
    // }

}