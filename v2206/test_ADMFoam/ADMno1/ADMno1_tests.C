// testing functions 

#include "ADMno1.H"

Foam::ADMno1::ADMno1
(
    volScalarField& T,
    const fvMesh& mesh,
    const IOdictionary& ADMno1Dict
)
:
    IOdictionary(ADMno1Dict),
    // DEBUG =======================================================
    Qin_
    (
        "Qin", 
        dimVolume/dimTime,
        // ADMno1Dict.lookupOrDefault("qin", 0.00)
        ADMno1Dict.lookupOrDefault("qin", 178.4674) // benchmark 
    ),
    Vgas_
    (
        "Vgas", 
        dimVolume,
        300
    ),
    Vliq_
    (
        "Vliq", 
        dimVolume, 
        3400
    ),
    Vfrac_test
    (
        IOobject
        (
            "Vfrac_test",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "Vfrac_testDefault", 
            dimless, 
            Zero
        )
    ),
    Vgas_test
    (
        IOobject
        (
            "Vgas_test",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "Vgas_testDefault", 
            dimVolume, 
            Zero
        )
    ),
    Ptotal_incell
    (
        IOobject
        (
            "Ptotal_incell",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
           "Ptotal_incellDefault", 
            dimPressure,
            Zero
        )
    ),
    rhoGas_test
    (
        IOobject
        (
            "rhoGas_test",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
           "rhoGas_testDefault", 
            dimDensity,
            Zero
        )
    ),
    vDotList_test(2),
    vDotGas_test
    (
        IOobject
        (
            "vDots_test",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            dimless/dimTime, 
            Zero
        )
    ),
    amplifier
    (
        ADMno1Dict.lookupOrDefault("amp", 1.0)
    ),
    // =============================================================
    para_
    (
        ADMno1Dict.get<word>("mode")
    ),
    Sc_
    (
        ADMno1Dict.lookupOrDefault("Sc", 1.0)
    ),
    Sct_
    (
        ADMno1Dict.lookupOrDefault("Sct", 0.2)
    ),
    R_
    (
        ADMno1Dict.lookupOrDefault("R", 0.083145 / para_.kTOK())
    ),
    KP_
    (
        ADMno1Dict.lookupOrDefault("Kpip", 5e4 / para_.BTOP())
    ),
    Vfrac_
    (
        ADMno1Dict.lookupOrDefault("Vfrac", 0.0882353) // 300/3400
    ), 
    Pvap_
    (
        "Pvap", 
        dimPressure,
        ADMno1Dict.lookupOrDefault("Pvap", para_.BTOP() * 0.0313)
    ),
    Pext_
    (
        "Pext", 
        dimPressure,
        ADMno1Dict.lookupOrDefault("Pext", para_.BTOP() * 1.013)
    ),
    Pgas_
    (
        IOobject
        (
            "Pgas",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        Pext_
    ),
    // TopDummy_(T), // DEBUG MULTI
    TopDummy_
    (
        IOobject
        (
            "TopDummy",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "TopDummy", 
            dimless, 
            308.15
        )
    ),
    fac_
    (
        IOobject
        (
            "fac",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "facDefault", 
            dimless, 
            ADMno1Dict.lookupOrDefault("fac", 1)
        )
    ),
    KHh2_(fac_),
    KHch4_(fac_),
    KHco2_(fac_),
    KaW_(fac_),
    KaIN_(fac_),                    
    Kaco2_(fac_),
    pH_
    (
        IOobject
        (
            "pH",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "pHDefault", 
            dimless, 
            ADMno1Dict.lookupOrDefault("pH", 7.26)
        )
    ),
    ShP_
    (
        IOobject
        (
            "ShP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
           "ShPDefault", 
            dimMass,
            para_.Pini()
        )
    ),
    Scat_
    (
        "Scat",
        dimMass/dimVolume, //TODO
        ADMno1Dict.lookupOrDefault("Scat", 0.00)
    ),
    San_
    (
        "San",
        dimMass/dimVolume, //TODO
        ADMno1Dict.lookupOrDefault("San", 0.0052 * para_.MTOm())
    ),
    tc_
    (
        "timeScale",
        dimTime, //TODO
        One
    )
{

    Info<< "\nSelecting ADM no1 operation mode " << ADMno1Dict.get<word>("mode") << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Main substances concentration initialization

    Info<< "Reading ADM no1 initial concentrations for soluables" << endl;

    label iNames = 0;
    label nSpecies = namesSoluable.size() + namesParticulate.size();

    YPtrs_.resize(nSpecies);

    forAll(namesSoluable, i)
    {
        YPtrs_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    namesSoluable[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }

    iNames += namesSoluable.size();

    //- Read particulates

    Info<< "Reading ADMno1 initial concentrations for particulates" << endl;

    forAll(namesParticulate, i)
    {
        YPtrs_.set
        (
            i + iNames,
            new volScalarField
            (
                IOobject
                (
                    namesParticulate[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }

    //- Initializing derivatives

    dYPtrs_.resize(nSpecies);

    for (label i = 0; i < nSpecies; i++)
    {
        dYPtrs_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    "d" + YPtrs_[i].name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    YPtrs_[0].dimensions()/dimTime, 
                    Zero
                )
            )
        );
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Gaseuoses initialization

    Info<< "Reading ADM no1 initial concentrations for gaseuoses" << endl;

    GPtrs_.resize(namesGaseous.size());

    forAll(namesGaseous, i)
    {
        GPtrs_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    namesGaseous[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    namesGaseous[i] + "Default", 
                    YPtrs_[0].dimensions(),
                    para_.Gini(i)
                )
            )
        );
    }
    
    //- Initializing derivatives

    dGPtrs_.resize(namesGaseous.size());

    for (label i = 0; i < namesGaseous.size(); i++)
    {
        dGPtrs_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    "d" + GPtrs_[i].name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    GPtrs_[0].dimensions()/dimTime, 
                    Zero
                )
            )
        );
    }


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Medians initialization

    Info<< "Initializing concentrations for medians" << endl;

    MPtrs_.resize(namesMedians.size());

    forAll(namesMedians, i)
    {
        MPtrs_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    namesMedians[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    namesMedians[i] + "Default", 
                    YPtrs_[0].dimensions(),
                    para_.Mini(i)
                )
            )
        );  
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Ions initialization

    IOPtrs_.resize(2);

    IOPtrs_.set
    (
        0,
        new volScalarField
        (
            IOobject
            (
                "Scat",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
               "Scat", 
                YPtrs_[0].dimensions(),
                ADMno1Dict.lookupOrDefault("Scat", 0.00)
            )
        )
    );

    IOPtrs_.set
    (
        1,
        new volScalarField
        (
            IOobject
            (
                "San",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),                             
            mesh,
            dimensionedScalar
            (
               "San", 
                YPtrs_[0].dimensions(),
                ADMno1Dict.lookupOrDefault("San", 0.0052 * para_.MTOm())
            )
        )
    );

    
    //- Initializing derivatives

    dIOPtrs_.resize(namesIons.size());

    for (label i = 0; i < namesIons.size(); i++)
    {
        dIOPtrs_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    "d" + IOPtrs_[i].name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    IOPtrs_[0].dimensions()/dimTime, 
                    Zero
                )
            )
        );
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //-  Medians initialization

    Info<< "Initializing concentrations for electrolytes" << endl;

    EPtrs_.resize(namesElectrolytes.size());

    forAll(namesElectrolytes, i)
    {
        EPtrs_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    namesElectrolytes[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ, // READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    YPtrs_[0].dimensions(),
                    para_.Eini(i)
                )
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Inhibition coeffs initialization

    IPtrs_.resize(8);

    for (int i = 0; i < 8; i++)
    {
        IPtrs_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Inh" + Foam::name(i),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE// TODO: choose if writing
                ),
                mesh,
                dimensionedScalar
                (
                    dimless, 
                    Zero
                )
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Kinetic rate initialization
    
    KRPtrs_.resize(19);

    for (int i = 0; i < 19; i++)
    {
        KRPtrs_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RRs" + Foam::name(i),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    dYPtrs_[0].dimensions(), 
                    Zero
                )
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Gas trasfer rate initialization
    
    GRPtrs_.resize(3);

    for (int i = 0; i < 3; i++)
    {
        GRPtrs_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    "GRs" + Foam::name(i),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    dGPtrs_[0].dimensions(), 
                    Zero
                )
            )
        );
    }

    // testing

    GPtrs_test.resize(namesGaseous.size());

    forAll(namesGaseous, i)
    {
        GPtrs_test.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    namesGaseous[i] + "_test",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    namesGaseous[i] + "_test_Default", 
                    YPtrs_[0].dimensions(),
                    Zero
                )
            )
        );
    }

    GRPtrs_test.resize(3);

    for (int i = 0; i < 3; i++)
    {
        GRPtrs_test.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    "GRs_test_" + Foam::name(i),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    dGPtrs_[0].dimensions(), 
                    Zero
                )
            )
        );
    }

    dGPtrs_test.resize(namesGaseous.size());

    for (label i = 0; i < namesGaseous.size(); i++)
    {
        dGPtrs_test.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    "d" + GPtrs_[i].name() + "_test",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    GPtrs_[0].dimensions()/dimTime, 
                    Zero
                )
            )
        );
    }

    vDotList_test.set
    (
        0,
        "gas",
        new volScalarField::Internal
        (
            IOobject
            (
                "vDotsG",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                dimless/dimTime, 
                Zero
            )
        )
    );

    vDotList_test.set
    (
        1,
        "sludge",
        new volScalarField::Internal
        (
            IOobject
            (
                "vDotsX",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                dimless/dimTime, 
                Zero
            )
        )
    );

    // vDotPtrs_test.resize(3);

    // for (int i = 0; i < 3; i++)
    // {
    //     vDotPtrs_test.set
    //     (
    //         i,
    //         new volScalarField::Internal
    //         (
    //             IOobject
    //             (
    //                 "vDots_test_" + Foam::name(i+1),
    //                 mesh.time().timeName(),
    //                 mesh,
    //                 IOobject::NO_READ,
    //                 IOobject::NO_WRITE
    //                 // IOobject::AUTO_WRITE
    //             ),
    //             mesh,
    //             dimensionedScalar
    //             (
    //                 dimVolume/dimTime, 
    //                 Zero
    //             )
    //         )
    //     );
    // }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // reset dimensions 
    para_.setParaDim(YPtrs_[0].dimensions());
    MPtrs_[0].ref() = YPtrs_[9] - EPtrs_[4]; // Sco2 = SIC - Shco3N

    ShP_.dimensions().reset(YPtrs_[0].dimensions());
    Scat_.dimensions().reset(YPtrs_[0].dimensions());
    San_.dimensions().reset(YPtrs_[0].dimensions());

    TopDummy_.dimensions().reset(dimless);

    KHh2_.dimensions().reset(para_.KH().h2.dimensions());
    KHch4_.dimensions().reset(para_.KH().ch4.dimensions());
    KHco2_.dimensions().reset(para_.KH().co2.dimensions());

    Kaco2_.dimensions().reset(para_.Ka().co2.dimensions());
    KaIN_.dimensions().reset(para_.Ka().IN.dimensions());
    KaW_.dimensions().reset(para_.Ka().W.dimensions());

    nIaa_ = 3.0 / (para_.pHL().ULaa - para_.pHL().LLaa);  // aa
    nIac_ = 3.0 / (para_.pHL().ULac - para_.pHL().LLac);  // ac
    nIh2_ = 3.0 / (para_.pHL().ULh2 - para_.pHL().LLh2);  // h2
    
    // DEBUG MULTI
    calcThermal(T);
}



void Foam::ADMno1::gasTest
(
    const volScalarField& T, // DEBUG MULTI
    const volScalarField &alphaLiq, 
    const volScalarField& Ptotal    
)
{
    volScalarField Ph2o_incell
    (
        IOobject
        (
            "Ph2o_incell",
            fac_.mesh().time().timeName(),
            fac_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fac_.mesh(),
        dimensionedScalar
        (
           "Ph2o_incellDefault", 
            GPtrs_[0].dimensions(),
            Zero
        )
    );

    volScalarField Pgas_incell
    (
        IOobject
        (
            "Pgas_incell",
            fac_.mesh().time().timeName(),
            fac_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fac_.mesh(),
        dimensionedScalar
        (
           "Pgas_incellDefault", 
            GPtrs_[0].dimensions(),
            Zero
        )
    );

    // TODO!!!: This needs to be reviewed since now we have multiphase and due to the volume fraction 
    //          of the gas phase, the concentration might differ from the origial calculation
    //          Consult solvers like [icoReactingMultiphaseInterFoam]
    // NOTE:    The concentration in the liquid phase doesn't not change with the phase vol fraction
    //          
    GRPtrs_test[0] = // <- in dimension of [kg * m-3 * s-1]
    (                // TODO: may need correction from saturation pressure
        para_.DTOS() * para_.kLa()                                 // GPtrs is gas fraction based concentration for now !!
      * (YPtrs_[7].internalField() - R_ * TopDummy_.internalField() * GPtrs_test[0].internalField() * KHh2_)
    );

    GRPtrs_test[1] = // <- in dimension of [kg * m-3 * s-1]
    (                // TODO: may need correction from saturation pressure
        para_.DTOS() * para_.kLa()                                 // GPtrs is gas fraction based concentration for now !!
      * (YPtrs_[8].internalField() - R_ * TopDummy_.internalField() * GPtrs_test[1].internalField() * KHch4_)
    );

    GRPtrs_test[2] = // <- in dimension of [mol * m-3 * s-1]
    (                // TODO: may need correction from saturation pressure
        para_.DTOS() * para_.kLa()                                 // GPtrs is gas fraction based concentration for now !!
      * (MPtrs_[0].internalField() - R_ * TopDummy_.internalField() * GPtrs_test[2].internalField() * KHco2_) 
    );// ^ Sco2 instead of SIC

    // GRPtrs_test[0] = // <- in dimension of [kg * m-3 * s-1]
    // (                // TODO: may need correction from saturation pressure
    //     para_.DTOS() * para_.kLa()
    //   * (YPtrs_[7].internalField() - Ph2 * KHh2_)
    // );

    // GRPtrs_test[1] = // <- in dimension of [kg * m-3 * s-1]
    // (                // TODO: may need correction from saturation pressure
    //     para_.DTOS() * para_.kLa()
    //   * (YPtrs_[8].internalField() - Pch4 * KHch4_)
    // );

    // GRPtrs_test[2] = // <- in dimension of [mol * m-3 * s-1]
    // (                // TODO: may need correction from saturation pressure
    //     para_.DTOS() * para_.kLa()                                    
    //   * (MPtrs_[0].internalField() - Pco2 * KHco2_) 
    // );// ^ Sco2 instead of SIC


    // ==================================================================================
    // field of cell volume for mesh 
    const dimensionedScalar TSat(dimTemperature, 310);
    const dimensionedScalar T0(dimTemperature, Zero);
    scalarField volMeshField = GPtrs_[0].mesh().V().field();

    // density of the bulk gas phase
    // rho = mi / V = Mi * ni / V = M * P / (R * T)
    // rhoGas_test.field() = 
    // (
    //     (Ptotal.field() / R_ * TopDummy_.internalField())
    //   * (
    //         (GPtrs_[0].internalField() / 8.0  + GPtrs_[1].internalField() / 4.0 + GPtrs_[2].internalField() * 0.044) 
    //       / (GPtrs_[0].internalField() / 16.0  + GPtrs_[1].internalField() / 64.0 + GPtrs_[2].internalField())
    //     )

    // vDotList_test["gas"].field() = // <-- check dimensions for Ptotal for multiphase
    // (  
    //     // amplifier * max(T - TSat, T0) // DEBUG MULTI (K) just to trigger phase change on the bottom
    //     amplifier * (T - TSat) / TSat * pos(T - TSat)
    // );
    // );



    // TODO: maybe just return a single vDot instead Gh2, Gch4 and Gco2 sepreately?
    // (m3) [Pa * K-1 * m3 * mol-1] * K * Pa-1 * [mol * m-3 * s-1]
    vDotList_test["gas"].field() = // <-- check dimensions for Ptotal for multiphase
    (  
        max
        (
            (/*volMeshField * */ R_ * TopDummy_.internalField() / Ptotal.field())                
          * (
                (para_.MTOm() * GRPtrs_test[2].field())
              + (para_.MTOm() * GRPtrs_test[0].field() / 16.0) // converting from kgCOD/m3 to mol/m3 (1kg H2 needs 8kg O2)                       
              + (para_.MTOm() * GRPtrs_test[1].field() / 64.0) // converting from kgCOD/m3 to mol/m3 (1kg CH4 needs 4kg O2)                  
            )
            ,
            scalar(0)
        ) * scalar(1e-5)
    );

    vDotGas_test.field() = vDotList_test["gas"].field();
    

    // ==================================================================================
    // bar
    // Ph2o_incell.field() = para_.KH().h2o * exp(5290.0 * fac_ * 100 * R_);
    // Pgas_incell.field() = (GPtrs_[0] / 16.0 + GPtrs_[1] / 64.0 + GPtrs_[2]) * R_ * TopDummy_;
    // Vfrac_test.field() = 100000*(Ph2o_incell.field() + Pgas_incell.field()) / Ptotal.field();

    // Vfrac_test = Pgas_incell.field();

    // Pgas_incell.field() = 
    // (
    //     (para_.MTOm() * GPtrs_[0] / 16.0 + para_.MTOm() * GPtrs_[1] / 64.0 + GPtrs_[2]) * R_ * TopDummy_ 
    //   + Ph2o_incell
    // );
    // Ph2o_incell.field() = Pvap_ * exp(5290.0 * fac_ * R_);
    // Ptotal_incell.field() = Ph2o_incell.field() + Pgas_incell.field();

    // Vfrac_test.field() = 100000 * Pgas_incell.field() / Ptotal.field();

    // volScalarField volMeshField = GPtrs_[0].mesh().V();            
    // volScalarField GRMass = GRPtrs_test * volMeshField;
    // volScalarField GRMolar = GRMass / molarMass;
    // volScalarField volGasRate = GRMolar * R * T / P; <- (P_rgh or P?)
    

    // ==================================================================================

    // // field of cell volume for mesh 
    // scalarField volMeshField = GPtrs_[0].mesh().V().field();            

    // scalarField volGas = volMeshField / (1.0 + (1.0/Vfrac_test));
    // scalarField volLiq = volMeshField / (1.0 + Vfrac_test);

    // Info<< "volGas: " << volGas << ", volLiq: " << volLiq << endl;

    // dGPtrs_test[0].field() = GRPtrs_test[0].field() * volLiq / volGas;
    // dGPtrs_test[1].field() = GRPtrs_test[1].field() * volLiq / volGas;
    // dGPtrs_test[2].field() = GRPtrs_test[2].field() * volLiq / volGas;

    // ==================================================================================

    // volScalarField moleRate_test_0 = (para_.MTOm() * GRPtrs_test[0].field() / 16.0) * volMeshField;
    // volScalarField moleRate_test_1 = (para_.MTOm() * GRPtrs_test[1].field() / 64.0) * volMeshField;
    // volScalarField moleRate_test_2 = GRPtrs_test[2].field() * volMeshField;

    // vDotPtrs_test[0].field() = moleRate_test_0.field() * R_ * TopDummy_.internalField() / Ptotal.field();
    // vDotPtrs_test[1].field() = moleRate_test_1.field() * R_ * TopDummy_.internalField() / Ptotal.field();
    // vDotPtrs_test[2].field() = moleRate_test_2.field() * R_ * TopDummy_.internalField() / Ptotal.field();
}




//- Functions for Sh2 calculations
volScalarField::Internal Foam::ADMno1::fSh2
(
    const surfaceScalarField &flux,
    const volScalarField &alphaLiq, 
    volScalarField &Sh2Temp
)
{
    volScalarField::Internal I_h2fa = calcInhibition // h2_fa
    (
        Sh2Temp,
        para_.KI().h2fa
    );

    volScalarField::Internal I_h2c4 = calcInhibition // h2_c4
    (
        Sh2Temp,
        para_.KI().h2c4
    );

    volScalarField::Internal I_h2pro = calcInhibition // h2_pro
    (
        Sh2Temp,
        para_.KI().h2pro
    );

    PtrList<volScalarField::Internal> KRPtrs_temp = KRPtrs_;

    KRPtrs_temp[6] = KRPtrs_temp[6]/IPtrs_[4]*I_h2fa;
    KRPtrs_temp[7] = KRPtrs_temp[7]/IPtrs_[5]*I_h2c4;
    KRPtrs_temp[8] = KRPtrs_temp[8]/IPtrs_[5]*I_h2c4;
    KRPtrs_temp[9] = KRPtrs_temp[9]/IPtrs_[6]*I_h2pro;

    KRPtrs_temp[11] = calcRho
    (
        para_.kDec().m_h2,
        Sh2Temp,
        para_.KS().h2,
        YPtrs_[22], // Xh2
        IPtrs_[2] * IPtrs_[3] // Iphh2*IIN
    );

    // volScalarField conv(fvc::div(flux, Sh2Temp));
    // volScalarField conv(fvc::div(flux, Sh2Temp)) / Vcell;
    volScalarField conv = para_.DTOS() * (Qin_/Vliq_) * (para_.INFLOW(7) - Sh2Temp);
    volScalarField::Internal GRSh2Temp = 
    (
        para_.DTOS() * para_.kLa() 
      * (Sh2Temp.internalField() - R_ * TopDummy_.internalField() * GPtrs_[0].internalField() * KHh2_)
    );

    //     reaction + convection - fGasRhoH2(paraPtr, Sh2);
    return concPerComponent(7, KRPtrs_temp) + conv - GRSh2Temp;
}


volScalarField::Internal Foam::ADMno1::dfSh2
(
    const surfaceScalarField &flux,
    const volScalarField &alphaLiq, 
    volScalarField &Sh2Temp
)
{
    volScalarField::Internal dI_h2fa = dCalcInhibition // h2_fa
    (
        Sh2Temp,
        para_.KI().h2fa
    );

    volScalarField::Internal dI_h2c4 = dCalcInhibition // h2_c4
    (
        Sh2Temp,
        para_.KI().h2c4
    );

    volScalarField::Internal dI_h2pro = dCalcInhibition // h2_pro
    (
        Sh2Temp,
        para_.KI().h2pro
    );

    PtrList<volScalarField::Internal> dKRPtrs_temp = KRPtrs_;
    forAll(dKRPtrs_temp, i)
    {
        dKRPtrs_temp[i].dimensions().reset(KRPtrs_[0].dimensions()/Sh2Temp.dimensions());
    }

    dKRPtrs_temp[6] = KRPtrs_[6]/IPtrs_[4]*dI_h2fa;
    dKRPtrs_temp[7] = KRPtrs_[7]/IPtrs_[5]*dI_h2c4;
    dKRPtrs_temp[8] = KRPtrs_[8]/IPtrs_[5]*dI_h2c4;
    dKRPtrs_temp[9] = KRPtrs_[9]/IPtrs_[6]*dI_h2pro;

    dKRPtrs_temp[11] = 
    (
        para_.kDec().m_h2 * YPtrs_[22].internalField() * IPtrs_[2] * IPtrs_[3] * para_.KS().h2
      / ((para_.KS().h2 + Sh2Temp.internalField()) * (para_.KS().h2 + Sh2Temp.internalField()))
    );

    // volScalarField dConv(fvc::div(flux));
    // volScalarField dConv(fvc::div(flux)) / Vcell;
    dimensionedScalar dConv = - para_.DTOS() * (Qin_/Vliq_);
    dimensionedScalar dGRSh2Temp = para_.DTOS() * para_.kLa();

    //     dReaction + dConvection - dfGasRhoH2(paraPtr, Sh2);
    return concPerComponent(7, dKRPtrs_temp) + dConv - dGRSh2Temp;
}


void Foam::ADMno1::calcSh2
(
    const surfaceScalarField &flux,
    const volScalarField &alphaLiq
)
{
    //TODO: IO dictionary for these parameters
    scalar tol = 1e-12;
    label nIter = 1e3;
    label i = 0;

    // initial value of x, E and dEdx
    volScalarField x = YPtrs_[7];   // x = Sh2
    volScalarField::Internal E = YPtrs_[7].internalField();   // E = dSh2/dt
    volScalarField::Internal dE = YPtrs_[7].internalField();  // dE = (dSh2/dt)/dSh2

    do
    {
        E.field() = fSh2(flux, x).field();
        dE.field() = dfSh2(flux, x).field();
        x.field() = x.field() - E.field()/dE.field();
        // false check
        // if( min(x.field()) < 0 )
        // {
        //     std::cerr << nl << "--> FOAM FATAL IO ERROR:" << nl
        //               << "Sh2 concentration below Zero\n";
        //     std::exit(1);
        //     break;
        // }
        // Info<< max(x.field()) << endl;
        i++;
    }
    while
    (
        max(mag(E.field())) > tol &&
        i < nIter
    );

    if( min(x.field()) < 0 )
    {
        x.field() = 0.0*x.field() + 1e-16;
    }

    Info<< "Newton-Raphson:\tSolving for Sh2" 
        << ", min Sh2: " << min(x.field()) 
        << ", max Sh2: " << max(x.field()) 
        << ", No Interations " << i << endl;

    // Sh2
    YPtrs_[7].ref() = x;
}


void Foam::ADMno1::dYUpdate
(
    const surfaceScalarField &flux,
    const volScalarField &alphaLiq
)
{
    for (label j = 0; j < 7; j++)
    {
        volScalarField::Internal inOutFlow(dYPtrs_[0]);
        inOutFlow = para_.DTOS() * (Qin_/Vliq_) * (para_.INFLOW(j) - YPtrs_[j]);
        dYPtrs_[j] = concPerComponent(j, KRPtrs_) + inOutFlow;
        // dYPtrs_[j] = concPerComponent(j, KRPtrs_);
    }

    for (label j = 8; j < YPtrs_.size(); j++)
    {
        volScalarField::Internal inOutFlow(dYPtrs_[0]);
        inOutFlow = para_.DTOS() * (Qin_/Vliq_) * (para_.INFLOW(j) - YPtrs_[j]);
        dYPtrs_[j] = concPerComponent(j, KRPtrs_) + inOutFlow;
        // dYPtrs_[j] = concPerComponent(j, KRPtrs_);
    }

    dIOPtrs_[0] = para_.DTOS() * (Qin_/Vliq_) * (para_.INFLOW(24) - IOPtrs_[0]);  // Scat
    dIOPtrs_[1] = para_.DTOS() * (Qin_/Vliq_) * (para_.INFLOW(25) - IOPtrs_[1]);  // San

    //- calculate with STOI and gas transer
    dYPtrs_[8] -= GRPtrs_[1]; // Sch4 - Gch4
    dYPtrs_[9] -= GRPtrs_[2]; // SIC - Gco2
}



void Foam::ADMno1::correct
(
    const surfaceScalarField &flux,
    const volScalarField &alphaLiq, 
    const volScalarField& T,
    const volScalarField& Ptotal
)
{
    //- Calculate thermal factor and adjust parameters
    // calcThermal(T); // DEBUG MULTI

    // testing <- not impacting the simulation for now
    // gasTest(T, Ptotal);
    gasTest(T, alphaLiq, Ptotal);

    //- Gas phase pressure
    gasPressure();

    //- Inhibition rates
    inhibitions();

    //- calculate raction rates
    kineticRate();

    //- calculate gas phase transfer rates
    gasPhaseRate();

    //- calculate dY with STOI
    // dYUpdate(flux, alphaLiq);
    dYUpdate(flux);

    //- calculate gas exit rates
    gasSourceRate();

    //- Acid-base calculations
    calcShp();

    //- Sh2 calculations
    calcSh2(flux);
}


Foam::PtrList<Foam::volScalarField>& Foam::ADMno1::G_test()
{
    return GPtrs_test;
}

Foam::PtrList<Foam::volScalarField::Internal>& Foam::ADMno1::dG_test()
{
    return dGPtrs_test;
}

tmp<fvScalarMatrix> Foam::ADMno1::RG_test
(
    label i
) const
{
    DimensionedField<scalar, volMesh> dG_test = dGPtrs_test[i];

        tmp<fvScalarMatrix> tSu
        (
            new fvScalarMatrix
            (
                GPtrs_test[i],
                dG_test.dimensions()*dimVolume
                // dimMass/dimTime // <- for compressible flow and uses fvm::ddt(rho, Yi)
            )
        );

    fvScalarMatrix& Su = tSu.ref();
    
    // https:// /documentation/guides/latest/api/fvMatrix_8C_source.html#l01708
    Su += dG_test; 

    return tSu;
};


