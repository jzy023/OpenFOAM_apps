


void Foam::multiphaseADMixture::solveAlphas
(
    const scalar cAlpha,
    const PtrListDictionary<volScalarField::Internal>& vDotList
)
{
    static label nSolves(-1);
    ++nSolves;

    const word alphaScheme("div(phi,alpha)");
    const word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic(mag(phi_/mesh_.magSf()));
    phic = min(cAlpha*phic, max(phic));

    // ------------------------------------------------------
    PtrList<surfaceScalarField> alphaPhiCorrs(phases_.size());
    int phasei = 0;

    for (phaseADM& alpha : phases_)
    {
        alphaPhiCorrs.set
        (
            phasei,
            new surfaceScalarField
            (
                "phi" + alpha.name() + "Corr",
                fvc::flux
                (
                    phi_,
                    alpha,
                    alphaScheme
                )
            )
        );

        surfaceScalarField& alphaPhiCorr = alphaPhiCorrs[phasei];

        for (phaseADM& alpha2 : phases_)
        {
            if (&alpha2 == &alpha) continue;

            surfaceScalarField phir(phic*nHatf(alpha, alpha2));

            alphaPhiCorr += fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha,
                alpharScheme
            );
        }

        MULES::limit
        (
            1.0/mesh_.time().deltaT().value(),
            geometricOneField(),
            alpha,
            phi_,
            alphaPhiCorr,
            zeroField(),
            zeroField(),
            oneField(),
            zeroField(),
            true
        );

        ++phasei;
    }

    MULES::limitSum(alphaPhiCorrs);

    rhoPhi_ = dimensionedScalar(dimMass/dimTime, Zero);

    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    // ------------------------------------------------------
    volScalarField::Internal alphaLiq = phases_.lookup("liquid")->internalField();

    // Set Su and Sp to zero
    for (phaseADM& alpha : phases_)
    {
        Su_[alpha.name()] = dimensionedScalar("Su", dimless/dimTime, Zero);
        Sp_[alpha.name()] = dimensionedScalar("Sp", dimless/dimTime, Zero);
    }

    // TODO: this works without addressing phase pairs because there are only 3 phases and
    //       there are only [gas -> liquid] and [liquid -> sludge] 
    for (phaseADM& alpha : phases_)
    {   
        // load vDot for gas and sludge generation
        if (alpha.name() == "liquid")
        {
            // TODO: to include Su['liquid'] for mass transfer to sludge 
            continue;
        }
        else
        {
            // adding Su of liquid to other phases
            Su_[alpha.name()] = vDotList[alpha.name()] * alphaLiq;
            // Sp_[alpha.name()] = [add implicity source term]

            // adding Su of other phases to liquid
            Su_["liquid"] -= vDotList[alpha.name()] * alphaLiq;
            // Sp_["liquid"] = [add implicity source term]
        }
    }

    // ------------------------------------------------------
    volScalarField divU
    (
        mesh_.moving()
      ? fvc::div(phi_ + mesh_.phi())
      : fvc::div(phi_)
    );

    phasei = 0;

    for (phaseADM& alpha : phases_)
    {
        surfaceScalarField& alphaPhi = alphaPhiCorrs[phasei];
        // alphaPhi += upwind<scalar>(mesh_, phi_).flux(alpha);

        const volScalarField::Internal& Su = Su_[alpha.name()];
        const volScalarField::Internal& Sp = Sp_[alpha.name()];

        // fvScalarMatrix alphaEqn
        // (
        //     fv::EulerDdtScheme<scalar>(mesh_).fvmDdt(alpha)
        //   + fv::gaussConvectionScheme<scalar>
        //     (
        //         mesh_,
        //         phi_,
        //         upwind<scalar>(mesh_, phi_)
        //     ).fvmDiv(phi_, alpha)
        //   - fvm::Sp(divU, alpha)
        //   ==
        //     Su 
        // //   + fvm::Sp(Sp, alpha)
        // );

        fvScalarMatrix alphaEqn
        (
            fv::EulerDdtScheme<scalar>(mesh_).fvmDdt(alpha)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh_,
                phi_,
                upwind<scalar>(mesh_, phi_)
            ).fvmDiv(phi_, alpha)
          ==
            Su
        );

        alphaEqn.boundaryManipulate(alpha.boundaryFieldRef());
        alphaEqn.solve();
        alphaPhi += alphaEqn.flux();

        MULES::explicitSolve
        (
            geometricOneField(),
            alpha,
            phi_,
            alphaPhi,
            Sp,
            Su,
            oneField(),
            zeroField()
        );

        rhoPhi_ += alphaPhi*alpha.rho();

        Info<< alpha.name() << " volume fraction, min, max = "
            << alpha.weightedAverage(mesh_.V()).value()
            << ' ' << min(alpha).value()
            << ' ' << max(alpha).value()
            << endl;

        sumAlpha += alpha;

        ++phasei;
    }

    Info<< "phaseADM-sum volume fraction, min, max = "
        << sumAlpha.weightedAverage(mesh_.V()).value()
        << ' ' << min(sumAlpha).value()
        << ' ' << max(sumAlpha).value()
        << endl;

    // Correct the sum of the phaseADM-fractions to avoid 'drift'
    volScalarField sumCorr(1.0 - sumAlpha);
    for (phaseADM& alpha : phases_)
    {
        alpha += alpha*sumCorr;
    }

    calcAlphas();
}

