Foam::tmp<Foam::fvScalarMatrix>  
Foam::multiphaseADMixture::volTransfer
(
    const volScalarField& p,
    const PtrListDictionary<volScalarField::Internal>& vDotList
)
{
    auto teqn = tmp<fvScalarMatrix>::New(p, dimVolume/dimTime);
    auto& eqn = teqn.ref();

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

    const volScalarField::Internal alphaLiq = phases_["liquid"].internalField();
    const volScalarField::Internal alphaLiqCrr
    (
        min(max(alphaLiq, scalar(0)), scalar(1))
    );

    for (phaseADM& alpha : phases_)
    {
        // load vDot for gas and sludge generation
        if (alpha.name() == "liquid")
        {
            // TODO: to include Su['liquid'] for mass transfer to sludge 
            continue;
        }

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

        if (phases_[alpha.name()].isRhoField())
        {    
            const volScalarField::Internal rhoField
            (
                rhoFieldsPtrs_[alpha.name()].internalField()
            );
        }

        //
        // Sp -=
        //     (vDotList["gas"] * alphaLiqCrr * phases_["liquid"].rho() * pos(alphaLiqCrr - scalar(0))) 
        //   * ((1 / rhoField) - (1 / phases_["liquid"].rho()));

        // Su -=
        //   - (vDotList["gas"] * alphaLiqCrr * phases_["liquid"].rho() * pos(alphaLiqCrr - scalar(0))) * 310
        //   * ((1 / rhoField) - (1 / phases_["liquid"].rho()));

        // TODO: need investigation
        Su +=
        (
            (vDotList[alpha.name()] * alphaLiqCrr * phases_["liquid"].rho() * pos(alphaLiqCrr - scalar(0)))
          * ((1 / rhoField) - (1 / phases_["liquid"].rho()))
        );
    }

    // eqn += fvm::Sp(Sp, p) + Su;
    eqn += Su;

    return teqn;
}