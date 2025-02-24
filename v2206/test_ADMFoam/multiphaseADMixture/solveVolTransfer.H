Foam::tmp<Foam::fvScalarMatrix>  
Foam::multiphaseADMixture::volTransfer
(
    const volScalarField& p,
    const PtrListDictionary<volScalarField::Internal>& vDotList
)
{
    auto teqn = tmp<fvScalarMatrix>::New(p, dimVolume/dimTime);
    auto& eqn = teqn.ref();

    // auto tSp = tmp<volScalarField>::New
    // (
    //     IOobject
    //     (
    //         "Sp",
    //         mesh_.time().timeName(),
    //         mesh_
    //     ),
    //     mesh_,
    //     dimensionedScalar(dimless/dimTime/dimPressure, Zero)
    // );
    // auto& Sp = tSp.ref();

    // auto tSu = tmp<volScalarField>::New
    // (
    //     IOobject
    //     (
    //         "Su",
    //         mesh_.time().timeName(),
    //         mesh_
    //     ),
    //     mesh_,
    //     dimensionedScalar(dimless/dimTime, Zero)
    // );
    // auto& Su = tSu.ref();


    auto tSp = tmp<volScalarField::Internal>::New
    (
        IOobject
        (
            "Sp",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime/dimPressure, Zero)
    );
    auto& Sp = tSp.ref();

    auto tSu = tmp<volScalarField::Internal>::New
    (
        IOobject
        (
            "Su",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    );
    auto& Su = tSu.ref();

    //
    const volScalarField::Internal alphaLiq = phases_["liquid"].internalField();

    const volScalarField::Internal alphaLiqCrr
    (
        min(max(alphaLiq, scalar(0)), scalar(1))
    );

    Su +=
    (
        (vDotList["gas"] * alphaLiqCrr * 1000 * pos(alphaLiqCrr - scalar(0)))
      * ((1 / 1.3) - (1 / 1000))
    );

    //
    eqn += fvm::Sp(Sp, p) + Su;
    return teqn;


    // forAllConstIters(this->totalPhasePairs(), iter)
    // {
    //     const phasePair& pair = iter()();

    //     const phaseModel& phase1 = pair.phase1();
    //     const phaseModel& phase2 = pair.phase2();

    //     // liquid to gas
    //     const phasePairKey key12
    //     (
    //         phase1.name(),
    //         phase2.name(),
    //         true
    //     );

    //     if (massTransferModels_.found(key12))
    //     {
    //         autoPtr<interfaceCompositionModel>& interfacePtr =
    //             massTransferModels_[key12];

    //         tmp<volScalarField> KSp =
    //             interfacePtr->KSp(interfaceCompositionModel::P, p);

    //         if (KSp.valid())
    //         {
    //             Sp -=
    //                 KSp.ref()
    //                *(
    //                     - this->coeffs(phase1.name())
    //                     + this->coeffs(phase2.name())
    //                 );
    //         }

    //         tmp<volScalarField> KSu =
    //             interfacePtr->KSu(interfaceCompositionModel::P, p);

    //         if (KSu.valid())
    //         {
    //             Su -=
    //                 KSu.ref()
    //                *(
    //                     - this->coeffs(phase1.name())
    //                     + this->coeffs(phase2.name())
    //                 );
    //         }

    //         // If linearization is not provided used full explicit
    //         if (!KSp.valid() && !KSu.valid())
    //         {
    //             Su -=
    //                 *dmdt_[key12]
    //                 *(
    //                     - this->coeffs(phase1.name())
    //                     + this->coeffs(phase2.name())
    //                 );
    //         }
    //     }

    //     // gas to liquid
    //     const phasePairKey key21
    //     (
    //         phase2.name(),
    //         phase1.name(),
    //         true
    //     );

    //     if (massTransferModels_.found(key21))
    //     {
    //         autoPtr<interfaceCompositionModel>& interfacePtr =
    //             massTransferModels_[key21];

    //         tmp<volScalarField> KSp =
    //             interfacePtr->KSp(interfaceCompositionModel::P, p);

    //         if (KSp.valid())
    //         {
    //             Sp +=
    //                 KSp.ref()
    //                *(
    //                     - this->coeffs(phase1.name())
    //                     + this->coeffs(phase2.name())
    //                 );
    //         }

    //         tmp<volScalarField> KSu =
    //             interfacePtr->KSu(interfaceCompositionModel::P, p);

    //         if (KSu.valid())
    //         {
    //             Su +=
    //                 KSu.ref()
    //                *(
    //                     - this->coeffs(phase1.name())
    //                     + this->coeffs(phase2.name())
    //                 );
    //         }

    //         // If linearization is not provided used full explicit
    //         if (!KSp.valid() && !KSu.valid())
    //         {
    //             Su +=
    //                 *dmdt_[key21]
    //                 *(
    //                     - this->coeffs(phase1.name())
    //                     + this->coeffs(phase2.name())
    //                 );
    //         }
    //     }

    // }

    // eqn += fvm::Sp(Sp, p) + Su;
    // return teqn;
}