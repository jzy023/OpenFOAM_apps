{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");
    
    // ------------------------------------------------------
    // TODO: declare as class member and init at constructor
    volScalarField::Internal Su1_
    (
        IOobject
        (
            "SuLiq",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimless/dimTime, Zero)
    );

    volScalarField::Internal Sp1_
    (
        IOobject
        (
            "SpLiq",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimless/dimTime, Zero)
    );

    volScalarField::Internal Su2_
    (
        IOobject
        (
            "SuGas",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimless/dimTime, Zero)
    );

    volScalarField::Internal Sp2_
    (
        IOobject
        (
            "SpGas",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimless/dimTime, Zero)
    );

    // ------------------------------------------------------
    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        interface.correct();
        surfaceScalarField phir("phir", phic*interface.nHatf());
        PtrList<surfaceScalarField> alphaPhiCorrs(2);

        alphaPhiCorrs.set
        (
            0,
            new surfaceScalarField
            (
                "phi" + alpha1.name() + "Corr",
                fvc::flux
                (
                    phi,
                    alpha1,
                    alphaScheme
                )
            )
        );
        
        alphaPhiCorrs.set
        (
            1,
            new surfaceScalarField
            (
                "phi" + alpha2.name() + "Corr",
                fvc::flux
                (
                    phi,
                    alpha2,
                    alphaScheme
                )
            )
        );

        surfaceScalarField& alphaPhiCorrs1 = alphaPhiCorrs[0];
        surfaceScalarField& alphaPhiCorrs2 = alphaPhiCorrs[1];

        // // tmp<surfaceScalarField> talphaPhiCorr1to2
        alphaPhiCorrs1 += fvc::flux
        (
            -fvc::flux(-phir, alpha2, alpharScheme),
            alpha1,
            alpharScheme
        );

        // // tmp<surfaceScalarField> talphaPhiCorr2to1
        alphaPhiCorrs2 += fvc::flux
        (
            -fvc::flux(-phir, alpha1, alpharScheme),
            alpha2,
            alpharScheme
        );

        // ------------------------------------------------------
        // Su1_ = dimensionedScalar("SuZero", dimless/dimTime, Zero);
        // Sp1_ = dimensionedScalar("SpZero", dimless/dimTime, Zero);
        // calculateSuSp();

        Su2_ = divU*alpha1;
        Sp2_ = mixture->vDotAlphal();

        // ------------------------------------------------------
        // const volScalarField::Internal& Su = Su1_;
        // const volScalarField::Internal& Sp = Sp1_;

        MULES::limit
        (
            1.0/mesh.time().deltaT().value(),
            geometricOneField(),
            alpha1,
            phi,
            alphaPhiCorrs[0],
            Sp1_,
            Su1_,
            oneField(),
            zeroField(),
            true
        );

        MULES::limit
        (
            1.0/mesh.time().deltaT().value(),
            geometricOneField(),
            alpha2,
            phi,
            alphaPhiCorrs[1],
            Sp2_,
            Su2_,
            oneField(),
            zeroField(),
            true
        );

        MULES::limitSum(alphaPhiCorrs);

        volScalarField sumAlpha
        (
            IOobject
            (
                "sumAlpha",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(dimless, Zero)
        );

        // ------------------------------------------------------
        // const volScalarField::Internal& Su = Su1_;
        // const volScalarField::Internal& Sp = Sp1_;

        // tmp<surfaceScalarField> talphaPhi;
        // surfaceScalarField& alphaPhi = alphaPhiCorrs[0];
        fvScalarMatrix alpha1Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha1)
        //   ==
        //     Su1_ + fvm::Sp(Sp1_, alpha1)
        );

        // fvScalarMatrix alpha2Eqn
        // (
        //     fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha2)
        //   + fv::gaussConvectionScheme<scalar>
        //     (
        //         mesh,
        //         phi,
        //         upwind<scalar>(mesh, phi)
        //     ).fvmDiv(phi, alpha2)
        // //   ==
        // //     Su2_ + fvm::Sp(Sp2_, alpha2)
        // );

        // alpha1Eqn.boundaryManipulate(alpha1.boundaryFieldRef());
        // alpha2Eqn.boundaryManipulate(alpha2.boundaryFieldRef());
        alpha1Eqn.solve();
        // alpha2Eqn.solve();

        // // talphaPhi = alpha1Eqn.flux();
        // alphaPhiCorrs[0] += alpha1Eqn.flux();
        // alphaPhiCorrs[1] += alpha1Eqn.flux();

        // MULES::explicitSolve
        // (
        //     geometricOneField(),
        //     alpha1,
        //     phi,
        //     alphaPhiCorrs[0],
        //     // talphaPhi.ref(),
        //     Sp1_,
        //     Su1_,
        //     oneField(),
        //     zeroField()
        // );

        // MULES::explicitSolve
        // (
        //     geometricOneField(),
        //     alpha2,
        //     phi,
        //     alphaPhiCorrs[1],
        //     Sp2_,
        //     Su2_,
        //     oneField(),
        //     zeroField()
        // );

        // if (aCorr == nAlphaCorr - 1)
        // {
        //     volScalarField sumAlpha
        //     (
        //         IOobject
        //         (
        //             "sumAlpha",
        //             mesh.time().timeName(),
        //             mesh
        //         ),
        //         mesh,
        //         dimensionedScalar(dimless, Zero)
        //     );

        //     // Reset rhoPhi
        //     rhoPhi = dimensionedScalar("rhoPhi", dimMass/dimTime, Zero);

        //     sumAlpha += alpha1;
        //     sumAlpha += alpha2;

        //     // Update rhoPhi
        //     rhoPhi += rho1 * alphaPhiCorrs[0];
        //     rhoPhi += rho2 * alphaPhiCorrs[1];

        //     Info<< "Phase-sum volume fraction, min, max = "
        //         << sumAlpha.weightedAverage(mesh.V()).value()
        //         << ' ' << min(sumAlpha).value()
        //         << ' ' << max(sumAlpha).value()
        //         << endl;

        //     volScalarField sumCorr(1.0 - sumAlpha);

        //     // liquid phase
        //     alpha1 += alpha1*sumCorr;

        //     Info<< alpha1.name() << " volume fraction = "
        //         << alpha1.weightedAverage(mesh.V()).value()
        //         << "  Min(alpha) = " << min(alpha1).value()
        //         << "  Max(alpha) = " << max(alpha1).value()
        //         << endl;

        //     // gas phase
        //     alpha2 += alpha2*sumCorr;

        //     Info<< alpha2.name() << " volume fraction = "
        //         << alpha2.weightedAverage(mesh.V()).value()
        //         << "  Min(alpha) = " << min(alpha2).value()
        //         << "  Max(alpha) = " << max(alpha2).value()
        //         << endl;
        // }

        // // interface.correct();

        // // // volScalarField alpha10("alpha10", alpha1);

        // // rhoPhi = talphaPhi()*(rho1 - rho2) + phi*rho2;

        // // // // Cache alphaPhi
        // // // alphaPhi10 = talphaPhi();

        // // // testing
        // // // limit alphas
        // // mixture->limitAlpha();
    }

    // ------------------------------------------------------

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
