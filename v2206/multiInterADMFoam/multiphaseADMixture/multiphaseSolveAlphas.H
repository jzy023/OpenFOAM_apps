void Foam::multiphaseADMixture::solveAlphas
(
    const scalar cAlpha
//    ,const PtrListDictionary<volScalarField::Internal>& vDotList
)
{
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

                // initialize rho as field
                const volScalarField::Internal rhoField
                (
                    IOobject
                    (
                        "rhoField",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    phases_[alpha.name()].rho()
                );
                
                // adding Su of liquid to other phases
                // Su_[alpha.name()] += 
                // (
                //     (vDotList[alpha.name()] * alphaLiqCrr * phases_["liquid"].rho() * pos(alphaLiqCrr - scalar(0))) 
                //   / rhoField
                // );

                // Sp_[alpha.name()] -=
                // (
                //     (vDotList[alpha.name()] * alphaLiqCrr * phases_["liquid"].rho() * pos(alphaLiqCrr - scalar(0)))
                //   * ((1 / rhoField) - (1 / phases_["liquid"].rho()))
                // );

                // // adding Su of other phases to liquid
                // Su_["liquid"] += 
                // (
                //   - (vDotList[alpha.name()] * alphaLiqCrr * phases_["liquid"].rho() * pos(alphaLiqCrr - scalar(0)))
                //   / phases_["liquid"].rho()
                // );

                // Sp_["liquid"] -=
                // (
                //   - (vDotList[alpha.name()] * alphaLiqCrr * phases_["liquid"].rho() * pos(alphaLiqCrr - scalar(0)))
                //   * ((1 / phases_["liquid"].rho()) - (1 / rhoField))
                // );
                
                const volScalarField::Internal alphaLiqLimited 
                (
                    max(min(alphaLiq, scalar(1)), scalar(0))
                );

                const volScalarField::Internal alphaGasLimited 
                (
                    max(min(alphaGas, scalar(1)), scalar(0))
                );

                const volScalarField::Internal coeffLiq = 
                (
                    (1 / phases_["liquid"].rho()) - (alphaLiqLimited * ((1 / phases_["liquid"].rho()) - (1 / rhoField)))
                );

                const volScalarField::Internal coeffGas = 
                (
                    (1 / rhoField) - (alphaGasLimited * ((1 / rhoField) - (1 / phases_["liquid"].rho())))
                );

                // adding Su of liquid to other phases
                // Sp_[alpha.name()] = [add implicity source term] !!!!
                Su_[alpha.name()] +=
                (
                    // vDotList[alpha.name()] 
                    mDotTest_ * (1/alphaW_) * actPatchCells_
                  * alphaLiqCrr * phases_["liquid"].rho() * pos(alphaLiqCrr - scalar(0))
                  * coeffGas
                );

                // adding Su of other phases to liquid
                // Sp_["liquid"] = [add implicity source term] !!!!
                // Su_["liquid"] -=
                // (
                //     // vDotList[alpha.name()] 
                //     mDotTest_ * (1/alphaW_) * actPatchCells_
                //   * alphaLiqCrr * phases_["liquid"].rho() * pos(alphaLiqCrr - scalar(0))
                //   * coeffLiq
                // );
            }
        }
    
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
                Sp,
                Su,
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
        phasei = 0;
        for (phaseADM& alpha : phases_)
        {

            const volScalarField::Internal& Su = Su_[alpha.name()];
            const volScalarField::Internal& Sp = Sp_[alpha.name()];
            surfaceScalarField& alphaPhi = alphaPhiCorrs[phasei];

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
                // rhoPhi_ += fvc::interpolate(phase.rho()) * phase.alphaPhi();
                rhoPhi_ += alpha.rho() * alpha.alphaPhi();
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

