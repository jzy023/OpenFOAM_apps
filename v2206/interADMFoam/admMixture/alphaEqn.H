{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phir("phir", phic*interface.nHatf());

    tmp<surfaceScalarField> talphaPhi;

    const volScalarField& vDotAlphal = mixture->vDotAlphal();

    if (MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha1)
          - fvm::Sp(divU, alpha1)
          ==
            fvm::Sp(vDotAlphal, alpha1)
        );

        alpha1Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        talphaPhi = alpha1Eqn.flux();
    }

    volScalarField alpha10("alpha10", alpha1);

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        tmp<surfaceScalarField> talphaPhiCorr
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            talphaPhiCorr.ref() -= talphaPhi();

            volScalarField alpha100("alpha100", alpha10);
            alpha10 = alpha1;

            MULES::correct
            (
                geometricOneField(),
                alpha1,
                talphaPhi(),
                talphaPhiCorr.ref(),
                vDotAlphal,
                (
                    divU*(alpha10 - alpha100)
                  - vDotAlphal*alpha10
                )(),
                oneField(),
                zeroField()
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                talphaPhi.ref() += talphaPhiCorr();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                talphaPhi.ref() += 0.5*talphaPhiCorr();
            }
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi,
                talphaPhiCorr.ref(),
                vDotAlphal,
                (divU*alpha1)(),
                oneField(),
                zeroField()
            );

            talphaPhi = talphaPhiCorr;
        }

        alpha2 = 1.0 - alpha1;
    }

    rhoPhi = talphaPhi()*(rho1 - rho2) + phi*rho2;

    // Cache alphaPhi
    alphaPhi10 = talphaPhi();

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
