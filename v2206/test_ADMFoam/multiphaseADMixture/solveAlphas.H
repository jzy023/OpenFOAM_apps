
void Foam::multiphaseADMixture::solveAlphas
(
    const scalar cAlpha,
    const PtrListDictionary<volScalarField::Internal>& vDotList
)
{
    // static label nSolves(-1);
    // ++nSolves;

    const dictionary& alphaControls = mesh_.solverDict("alpha");
    label nAlphaCorr(alphaControls.get<label>("nAlphaCorr"));

    const word alphaScheme("div(phi,alpha)");
    const word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic(mag(phi_/mesh_.magSf()));

    // ------------------------------------------------------
    PtrList<surfaceScalarField> alphaPhiCorrs(phases_.size());
    
    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
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

                phic = min(cAlpha*phic, max(phic));

                surfaceScalarField phir(phic*nHatf(alpha, alpha2));

                alphaPhiCorr += fvc::flux
                (
                    -fvc::flux(-phir, alpha2, alpharScheme),
                    alpha,
                    alpharScheme
                );
            }

            // Ensure that the flux at inflow BCs is preserved
            // TODO:

            ++phasei;
        }                                     

    // ------------------------------------------------------
        const volScalarField::Internal alphaLiq = phases_["liquid"].internalField();
        const volScalarField::Internal alphaGas = phases_["gas"].internalField();

        // Set Su and Sp to zero
        for (phaseADM& alpha : phases_)
        {
            Su_[alpha.name()] = dimensionedScalar("Su", dimless/dimTime, Zero);
            Sp_[alpha.name()] = dimensionedScalar("Sp", dimless/dimTime, Zero);
        }

        // Including divU 
        Su_["liquid"] +=
            fvc::div(phi_)*min(max(phases_["liquid"], scalar(0)), scalar(1));

        Su_["gas"] +=
            fvc::div(phi_)*min(max(phases_["gas"], scalar(0)), scalar(1));


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
                const volScalarField::Internal alphaLiqCrr
                (
                    min(max(alphaLiq, scalar(0)), scalar(1))
                );

                // adding Su of liquid to other phases
                Su_[alpha.name()] += 
                (
                    (vDotList[alpha.name()] * alphaLiqCrr * 1000 * pos(alphaLiqCrr - scalar(0))) 
                  / 1.3
                );

                Sp_[alpha.name()] -=
                (
                    (vDotList[alpha.name()] * alphaLiqCrr * 1000 * pos(alphaLiqCrr - scalar(0)))
                  * ((1 / 1.3) - (1 / 1000))
                );

                // adding Su of other phases to liquid
                Su_["liquid"] += 
                (
                  - (vDotList[alpha.name()] * alphaLiqCrr * 1000 * pos(alphaLiqCrr - scalar(0)))
                  / 1000
                );

                Sp_["liquid"] -=
                (
                  - (vDotList[alpha.name()] * alphaLiqCrr * 1000 * pos(alphaLiqCrr - scalar(0)))
                  * ((1 / 1000) - (1 / 1.3))
                );

                // const volScalarField::Internal alphaLiqLimited 
                // (
                //     max(min(alphaLiq, scalar(1)), scalar(0))
                // );

                // const volScalarField::Internal alphaGasLimited 
                // (
                //     max(min(alphaGas, scalar(1)), scalar(0))
                // );

                // const volScalarField::Internal coeffLiq = 
                // (
                //     (1 / 1000) - (alphaLiqLimited * ((1 / 1000) - (1 / 1.3)))
                // );

                // const volScalarField::Internal coeffGas = 
                // (
                //     (1 / 1.3) - (alphaGasLimited * ((1 / 1.3) - (1 / 1000)))
                // );

                // // adding Su of liquid to other phases
                // // Sp_[alpha.name()] = [add implicity source term] !!!!
                // Su_[alpha.name()] +=
                // (
                //     (vDotList[alpha.name()] * alphaLiqCrr * 1000 * pos(alphaLiqCrr - scalar(0)))
                //   * coeffGas
                // );

                // // adding Su of other phases to liquid
                // // Sp_["liquid"] = [add implicity source term] !!!!
                // Su_["liquid"] -=
                // (
                //     (vDotList[alpha.name()] * alphaLiqCrr * 1000 * pos(alphaLiqCrr - scalar(0)))
                //   * coeffLiq
                // );
            }

           // forAll(Su_, celli)
           // {
           //     if (Su_["gas"][celli] > 1e-16)
           //     {
           //         Info<< ">>> DEBUG!!!"
           //             << "\n>>> celli: " << celli
           //             << "\n>>> (Su[gas] - Su[liquid]): " << (Su_["gas"][celli] + Su_["liquid"][celli])
           //             << "\n>>> (coeffGas - coeffLiq): "  << (coeffGas[celli] - coeffLiq[celli])
           //             << "\n>>> coeffGas: "               << coeffGas[celli]
           //             << endl;

           //         break;
           //     }
           // }
        }

        // const volScalarField::Internal alphaLiqCrr
        // (
        //     min(max(alphaLiq, scalar(0)), scalar(1))
        // );

        // // adding Su of liquid to other phases
        // // Sp_[alpha.name()] = [add implicity source term] !!!!
        // Su_["gas"] += 
        // (
        //     (vDotList["gas"] * alphaLiqCrr * 1000 * pos(alphaLiqCrr - scalar(0))) 
        //   / 1.3
        // );

        // Sp_["gas"] -=
        // (
        //     (vDotList["gas"] * alphaLiqCrr * 1000 * pos(alphaLiqCrr - scalar(0)))
        //   * ((1 / 1.3) - (1 / 1000))
        // );

        // // adding Su of other phases to liquid
        // // Sp_["liquid"] = [add implicity source term] !!!!
        // Su_["liquid"] += 
        // (
        //   - (vDotList["gas"] * alphaLiqCrr * 1000 * pos(alphaLiqCrr - scalar(0)))
        //   / 1000
        // );

        // Sp_["liquid"] -=
        // (
        //   - (vDotList["gas"] * alphaLiqCrr * 1000 * pos(alphaLiqCrr - scalar(0)))
        //   * ((1 / 1000) - (1 / 1.3))
        // );

        // const volScalarField::Internal alphaLiqLimited 
        // (
        //     max(min(alphaLiq, scalar(1)), scalar(0))
        // );

        // const volScalarField::Internal alphaGasLimited 
        // (
        //     max(min(alphaGas, scalar(1)), scalar(0))
        // );

        // const volScalarField::Internal coeffLiq = 
        // (
        //     (1 / 1000) - alphaLiqLimited 
        //   * ((1 / 1000) - (1 / 1.3))
        // );

        // const volScalarField::Internal coeffGas = 
        // (
        //     (1 / 1.3) - alphaGasLimited 
        //   * ((1 / 1.3) - (1 / 1000))
        // );

        // // adding Su of liquid to other phases
        // // Sp_["gas"] = [add implicity source term] !!!!
        // Su_["gas"] +=
        // (
        //     (vDotList["gas"] * alphaLiqCrr * 1000 * pos(alphaLiqCrr - scalar(0)))
        //   * coeffGas
        // );

        // // adding Su of other phases to liquid
        // // Sp_["liquid"] = [add implicity source term] !!!!
        // Su_["liquid"] +=
        // (
        //   - (vDotList["gas"] * alphaLiqCrr * 1000 * pos(alphaLiqCrr - scalar(0)))
        //   * coeffLiq
        // );
    
    // ------------------------------------------------------
        // Limit phiAlphaCorr on each phase
        phasei = 0;
        for (phaseADM& alpha : phases_)
        {
            surfaceScalarField& alphaPhiCorr = alphaPhiCorrs[phasei];

            const volScalarField::Internal& Su = Su_[alpha.name()];
            const volScalarField::Internal& Sp = Sp_[alpha.name()];

            MULES::limit
            (
                1.0/mesh_.time().deltaT().value(),
                geometricOneField(),
                alpha,
                phi_,
                alphaPhiCorr,
                Sp, // zeroField(),
                Su, // zeroField(),
                oneField(),
                zeroField(),
                true
            );

            ++phasei;
        }

        MULES::limitSum(alphaPhiCorrs);

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
        // volScalarField divU
        // (
        //     mesh_.moving()
        // ? fvc::div(phi_ + mesh_.phi())
        // : fvc::div(phi_)
        // );

        phasei = 0;
        for (phaseADM& alpha : phases_)
        {

            const volScalarField::Internal& Su = Su_[alpha.name()];
            const volScalarField::Internal& Sp = Sp_[alpha.name()];
            surfaceScalarField& alphaPhi = alphaPhiCorrs[phasei];

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
            //   + fvm::Sp(Sp, alpha)
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
            //   + fvm::Sp(Sp, alpha)
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

            alpha.alphaPhi() = alphaPhi;

            ++phasei;
        }

        // Info<< ">>> DEBUG!!!"
        //     << "\n>>> alpha sum: " << max(phases_["gas"] + phases_["liquid"])
        //     << endl;

        // forAll(phases_["gas"], celli)
        // {
        //     if ((phases_["gas"][celli] + phases_["liquid"][celli]) > 1.1)
        //     {
        //         Info<< ">>> DEBUG!!!"
        //             << "\n>>> celli: " << celli
        //             << "\n>>> alpha sum: " << (phases_["gas"][celli] + phases_["liquid"][celli])
        //             // << "\n>>> (Su[gas] - Su[liquid]): " << (Su_["gas"][celli] + Su_["liquid"][celli])
        //             // << "\n>>> (coeffGas - coeffLiq): "  << (coeffGas[celli] - coeffLiq[celli])
        //             // << "\n>>> coeffGas: "               << coeffGas[celli]
        //             << endl;
         
        //         break;
        //     }
        // }

        if (acorr == nAlphaCorr - 1)
        {
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

            // Reset rhoPhi
            rhoPhi_ = dimensionedScalar("rhoPhi", dimMass/dimTime, Zero);

            for (phaseADM& alpha : phases_)
            {
                sumAlpha += alpha;

                // Update rhoPhi
                rhoPhi_ += alpha.rho() * alpha.alphaPhi();
                // rhoPhi_ += fvc::interpolate(alpha.rho()) * alpha.alphaPhi();
            }

            Info<< "Phase-sum volume fraction, min, max = "
                << sumAlpha.weightedAverage(mesh_.V()).value()
                << ' ' << min(sumAlpha).value()
                << ' ' << max(sumAlpha).value()
                << endl;

            volScalarField sumCorr(1.0 - sumAlpha);

            for (phaseADM& alpha : phases_)
            {
                alpha += alpha*sumCorr;

                Info<< alpha.name() << " volume fraction = "
                    << alpha.weightedAverage(mesh_.V()).value()
                    << "  Min(alpha) = " << min(alpha).value()
                    << "  Max(alpha) = " << max(alpha).value()
                    << endl;
            }
        }
    }
}

