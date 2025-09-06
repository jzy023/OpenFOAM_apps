/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |  
     \\/     M anipulation  |
-------------------------------------------------------------------------------
 
License
    This file is part of OpenFOAM.
 
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
 
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
 
    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
 
\*---------------------------------------------------------------------------*/

#include "ADMno1.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
 
// namespace Foam
// {
//     defineTypeNameAndDebug(ADMno1, 0);
//     // defineRunTimeSelectionTable(ADMno1, fvMesh);
// }
 
const Foam::word Foam::ADMno1::propertiesName("admno1Properties");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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
        // 100 // <<< Rosen et al.
        300
    ),
    Vliq_
    (
        "Vliq", 
        dimVolume, 
        3400
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
        ADMno1Dict.lookupOrDefault
        (
            "R", 
            0.083145 / para_.kTOK()
            // default digit number from Rosen paper with their dimensions
        ) 
    ),
    kLa_
    (
        "kLa",
        dimless/dimTime,
        ADMno1Dict.lookupOrDefault
        (
            "kLa",
            200.0
        ) 
    ),
    KP_
    (
        ADMno1Dict.lookupOrDefault
        (
            "Kpip", 
            5e4 / para_.BTOP()
            // default digit number from Rosen paper with their dimensions
        ) 
    ),
    Vfrac_
    (
        // ADMno1Dict.lookupOrDefault("Vfrac", 0.0294118) // 100/3400
        ADMno1Dict.lookupOrDefault("Vfrac", 0.0882353) // 300/3400
    ), 
    Pvap_
    (
        "Pvap", 
        dimPressure,
        ADMno1Dict.lookupOrDefault
        (
            "Pvap", 
            para_.BTOP() * 0.0313
            // default digit number from Rosen paper with their dimensions
        ) 
    ),
    Pext_
    (
        "Pext", 
        dimPressure,
        ADMno1Dict.lookupOrDefault
        (
            "Pext", 
            para_.BTOP() * 1.013
            // default digit number from Rosen paper with their dimensions
        ) 
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
    TopDummy_(T),
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
                    namesSoluable[i], // IOobject::groupName(namesSoluable[i]),
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

    // DEBUG
    Vfrac_ = (Vgas_/Vliq_).value();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
 
Foam::autoPtr<Foam::ADMno1> Foam::ADMno1::New
(
    volScalarField& T,
    const fvMesh& mesh
)
{
    IOdictionary ADMno1Dict
    (
        IOobject
        (
            propertiesName,
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    // TODO: do it properly!!! with virtual destructors and constructor hash tables 
    // new keywaord is not gonna last!
    ADMno1* reactionPtr = new ADMno1(T, mesh, ADMno1Dict);
    return autoPtr<ADMno1>(reactionPtr);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ADMno1::calcThermal
(
    volScalarField& T
)
{
    TopDummy_.field() = T.field();

    fac_ = (1.0 / para_.Tbase().value() - 1.0 / TopDummy_) / R_;
    
    KHh2_ = para_.KH().h2 * exp(-4180.0 * fac_);
    KHch4_ = para_.KH().ch4 * exp(-14240.0 * fac_);
    KHco2_ = para_.KH().co2 * exp(-19410.0 * fac_);
    
    Kaco2_ = para_.Ka().co2 * exp(7646.0 * fac_);
    KaIN_ = para_.Ka().IN * exp(51965.0 * fac_);
    KaW_ = para_.Ka().W * exp(55900.0 * fac_);

    // Info<< "KHco2: " << max(KHco2_.field()) << endl;
    // Info<< "Kaco2: " << max(Kaco2_.field()) << endl;
    // Info<< "KaIN_: " << max(KaIN_.field()) << endl;
    // Info<< "KaW_: " << max(KaW_.field()) << endl;

}


//- Functions for gas transfer calculations
void Foam::ADMno1::gasPressure()
{
    //- gas pressure
    volScalarField Ph2o
    (
        IOobject
        (
            "Ph2o",
            fac_.mesh().time().timeName(),
            fac_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fac_.mesh(),
        dimensionedScalar
        (
           "Ph2oDefault", 
            GPtrs_[0].dimensions(),
            Zero
        )
    );

    Ph2o.field() = Pvap_ * exp(5290.0 * fac_ * R_);
    Pgas_.field() = 
    (
        Ph2o + R_ * TopDummy_
     * (para_.MTOm() * GPtrs_[0] / 16.0 + para_.MTOm() * GPtrs_[1] / 64.0 + GPtrs_[2])
    );    
}


void Foam::ADMno1::gasPhaseRate()
{
    GRPtrs_[0] = 
    (
        para_.DTOS() * kLa_
      * (YPtrs_[7].internalField() - R_ * TopDummy_.internalField() * GPtrs_[0].internalField() * KHh2_)
    );

    GRPtrs_[1] = 
    (
        para_.DTOS() * kLa_
      * (YPtrs_[8].internalField() - R_ * TopDummy_.internalField() * GPtrs_[1].internalField() * KHch4_)
    );

    GRPtrs_[2] = 
    (
        para_.DTOS() * kLa_// Sco2 instead of SIC
      * (MPtrs_[0].internalField() - R_ * TopDummy_.internalField() * GPtrs_[2].internalField() * KHco2_)
    );

    // test -------------------------------------------------------------
    dimensionedScalar Sh2Ave = YPtrs_[7].weightedAverage(YPtrs_[0].mesh().V());
    dimensionedScalar Sch4Ave = YPtrs_[8].weightedAverage(YPtrs_[0].mesh().V());
    dimensionedScalar Sco2Ave = MPtrs_[0].weightedAverage(YPtrs_[0].mesh().V());

    Info<< ">>> Sh2Ave = " << Sh2Ave.value() << endl;
    Info<< ">>> Sch4Ave = " << Sch4Ave.value() << endl;
    Info<< ">>> Sco2Ave = " << Sco2Ave.value() << endl;

    forAll(GPtrs_, i)
    {
        dimensionedScalar GiAve = GPtrs_[i].weightedAverage(GPtrs_[i].mesh().V());
        dimensionedScalar GRiAve = GRPtrs_[i].weightedAverage(GPtrs_[i].mesh().V());

        Info<< ">>> "                 << GPtrs_[i].name()
            << " concentration = "    << GiAve.value()
            << ", generation rate = " << GRiAve.value() << endl;
    }
    // ------------------------------------------------------------------
}


void Foam::ADMno1::gasSourceRate()
{
    // field of cell volume for mesh 
    scalarField volMeshField = GPtrs_[0].mesh().V().field();            

    // particle scaled gas volume
    scalarField volGas = volMeshField / (1.0 + (1.0/Vfrac_));

    // particle scaled liquid volume
    scalarField volLiq = volMeshField / (1.0 + Vfrac_);

    // particle scaled pipe resistance
    volScalarField kp
    (
        IOobject
        (
            "kp",
            fac_.mesh().time().timeName(),
            fac_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fac_.mesh(),
        dimensionedScalar
        (
           "kp_Default", 
            dimless,
            Zero
        )
    );

    // TODO: might lead to error if Gas dimension is different
    kp.dimensions().reset(dimVolume/dimTime/dimPressure);

    // TODO: actual volume would have effect on normalized kp
    // (Vgas_ + Vliq_) <---- this would be calculated
    kp.field() = para_.DTOS() * KP_ * (volMeshField / (Vgas_ + Vliq_).value());

    //  volScalarField qGasLocal = kp * (Pgas - Pext_) * (Pgas / Pext_); 
    volScalarField qGasLocal = kp * (Pgas_ - Pext_);
    forAll( qGasLocal.field(), i )
    {
        if ( qGasLocal.field()[i] < 0.0 ) { qGasLocal.field()[i] = 1e-16; }
    }

    dGPtrs_[0].field() = 
    (
        (GRPtrs_[0].field() * volLiq / volGas) 
      - (GPtrs_[0].field() * qGasLocal.field() / volGas)
    );

    dGPtrs_[1].field() = 
    (
        (GRPtrs_[1].field() * volLiq / volGas) 
      - (GPtrs_[1].field() * qGasLocal.field() / volGas)
    );

    dGPtrs_[2].field() = 
    (
        (GRPtrs_[2].field() * volLiq / volGas) 
      - (GPtrs_[2].field() * qGasLocal.field() / volGas)
    );
};


//- Functions for Sh2 calculations
volScalarField::Internal Foam::ADMno1::fSh2
(
    const surfaceScalarField &flux,
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
    volScalarField conv = para_.DTOS() * (Qin_/Vliq_) * (para_.INFLOW(7) - Sh2Temp);
    volScalarField::Internal GRSh2Temp = 
    (
        para_.DTOS() * kLa_ 
      * (Sh2Temp.internalField() - R_ * TopDummy_.internalField() * GPtrs_[0].internalField() * KHh2_)
    );

    //     reaction + convection - fGasRhoH2(paraPtr, Sh2);
    return concPerComponent(7, KRPtrs_temp) + conv - GRSh2Temp;
}


volScalarField::Internal Foam::ADMno1::dfSh2
(
    const surfaceScalarField &flux,
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
    dimensionedScalar dConv = - para_.DTOS() * (Qin_/Vliq_);
    dimensionedScalar dGRSh2Temp = para_.DTOS() * kLa_;

    //     dReaction + dConvection - dfGasRhoH2(paraPtr, Sh2);
    return concPerComponent(7, dKRPtrs_temp) + dConv - dGRSh2Temp;
}


void Foam::ADMno1::calcSh2
(
    const surfaceScalarField &flux
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
        gMax(mag(E.field())) > tol &&
        i < nIter
    );

    x.field() = min(max(x.field(), scalar(1e-16)), x.field());

    // if
    // (
    //     min(x.field()) < 0
    // )
    // {
    //     x.field() = 0.0*x.field() + 1e-16;
    // }

    Info<< "Newton-Raphson:\tSolving for Sh2" 
        << ", min Sh2: " << min(x.field()) 
        << ", max Sh2: " << max(x.field()) 
        << ", No Interations " << i << endl;

    // Sh2
    YPtrs_[7].ref() = x;
}


//- Functions for inhibition and kinetic rates calculations
void Foam::ADMno1::inhibitions()
{
    //- Inhibiitons

    IPtrs_[0] = calcInhibitionHP // pH_aa
    (
        ShP_,
        para_.pHL().ULaa, 
        para_.pHL().LLaa,
        nIaa_
    );

    IPtrs_[1] = calcInhibitionHP // pH_ac
    (
        ShP_,
        para_.pHL().ULac, 
        para_.pHL().LLac,
        nIac_
    );

    IPtrs_[2] = calcInhibitionHP // pH_h2
    (
        ShP_,
        para_.pHL().ULh2, 
        para_.pHL().LLh2,
        nIh2_
    );


    IPtrs_[3] = 1.0 / (1.0 + (para_.KS().IN / YPtrs_[10]));
    
    IPtrs_[4] = calcInhibition // h2_fa
    (
        YPtrs_[7], // Sh2
        para_.KI().h2fa
    );

    IPtrs_[5] = calcInhibition // h2_c4
    (
        YPtrs_[7], // Sh2
        para_.KI().h2c4
    );

    IPtrs_[6] = calcInhibition // h2_pro
    (
        YPtrs_[7], // Sh2
        para_.KI().h2pro
    );

    IPtrs_[7] = calcInhibition // nh3
    (
        MPtrs_[1], // Snh3
        para_.KI().nh3
    );
}


void Foam::ADMno1::kineticRate()
{
    //- Kinetic rates

    KRPtrs_[0] = calcRho
    (
        para_.kDec().dis,
        YPtrs_[12] // Xc
    );

    KRPtrs_[1] = calcRho
    (
        para_.kDec().hyd_ch,
        YPtrs_[13] // Xch
    );

    KRPtrs_[2] = calcRho
    (
        para_.kDec().hyd_pr,
        YPtrs_[14] // Xpr
    );

    KRPtrs_[3] = calcRho
    (
        para_.kDec().hyd_li,
        YPtrs_[15] // Xli
    );

    KRPtrs_[4] = calcRho
    (
        para_.kDec().m_su,
        YPtrs_[0], // Ssu
        para_.KS().su,
        YPtrs_[16], // Xsu
        IPtrs_[0] * IPtrs_[3] // Iphaa*IIN
    );

    KRPtrs_[5] = calcRho
    (
        para_.kDec().m_aa,
        YPtrs_[1], // Saa
        para_.KS().aa,
        YPtrs_[17], // Xaa
        IPtrs_[0] * IPtrs_[3] // Iphaa*IIN
    );

    KRPtrs_[6] = calcRho
    (
        para_.kDec().m_fa,
        YPtrs_[2], // Sfa
        para_.KS().fa,
        YPtrs_[18], // Xfa
        IPtrs_[0] * IPtrs_[3] * IPtrs_[4] //Iphaa*IIN*Ih2fa
    );

    KRPtrs_[7] = calcRho
    (
        para_.kDec().m_c4,
        YPtrs_[3], // Sva
        para_.KS().c4,
        YPtrs_[19], // Xc4
        YPtrs_[4],  // Sbu
        IPtrs_[0] * IPtrs_[3] * IPtrs_[5] //Iphaa*IIN*Ih2c4
    );

    KRPtrs_[8] = calcRho
    (
        para_.kDec().m_c4,
        YPtrs_[4], // Sbu
        para_.KS().c4,
        YPtrs_[19], // Xc4
        YPtrs_[3],  // Sva
        IPtrs_[0] * IPtrs_[3] * IPtrs_[5] //Iphaa*IIN*Ih2c4
    );

    KRPtrs_[9] = calcRho
    (
        para_.kDec().m_pro,
        YPtrs_[5], // Spro
        para_.KS().pro,
        YPtrs_[20], // Xpro
        IPtrs_[0] * IPtrs_[3] * IPtrs_[6]  //Iphaa*IIN*Ih2pro
    );

    KRPtrs_[10] = calcRho
    (
        para_.kDec().m_ac,
        YPtrs_[6], // Sac
        para_.KS().ac,
        YPtrs_[21], // Xac
        IPtrs_[1] * IPtrs_[3] * IPtrs_[7] // Iphac*IIN*Inh3
    );

	// >>> in Rosen et al. implementation, no intermediate used for S_h2
    KRPtrs_[11] = calcRho
    (
        para_.kDec().m_h2,
        YPtrs_[7], // Sh2
        para_.KS().h2,
        YPtrs_[22], // Xh2
        IPtrs_[2] * IPtrs_[3] // Iphh2*IIN
    );

    KRPtrs_[12] = calcRho
    (
        para_.kDec().dec_xsu,
        YPtrs_[16] // Xsu
    );

    KRPtrs_[13] = calcRho
    (
        para_.kDec().dec_xaa,
        YPtrs_[17] // Xaa
    );

    KRPtrs_[14] = calcRho
    (
        para_.kDec().dec_xfa,
        YPtrs_[18] // Xfa
    );

    KRPtrs_[15] = calcRho
    (
        para_.kDec().dec_xc4,
        YPtrs_[19] // Xc4
    );

    KRPtrs_[16] = calcRho
    (
        para_.kDec().dec_xpro,
        YPtrs_[20] // Xpro
    );

    KRPtrs_[17] = calcRho
    (
        para_.kDec().dec_xac,
        YPtrs_[21] // Xac
    );

    KRPtrs_[18] = calcRho
    (
        para_.kDec().dec_xh2,
        YPtrs_[22] // Xh2
    );
}


void Foam::ADMno1::dYUpdate
(
    const surfaceScalarField &flux
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


//- Functions for acid base reaction calculations
volScalarField::Internal Foam::ADMno1::fShp
(
    volScalarField::Internal& ShpTemp
)
{
    EPtrs_[0] = fSion
    (
        para_.Ka().va,
        YPtrs_[3].internalField(),
        ShpTemp
    );

    EPtrs_[1] = fSion
    (
        para_.Ka().bu,
        YPtrs_[4].internalField(),
        ShpTemp
    ); 

    EPtrs_[2] = fSion
    (
        para_.Ka().pro,
        YPtrs_[5].internalField(),
        ShpTemp
    ); 

    EPtrs_[3] = fSion
    (
        para_.Ka().ac,
        YPtrs_[6].internalField(),
        ShpTemp
    );

    EPtrs_[4] = fSion
    (
        Kaco2_,
        YPtrs_[9].internalField(), // SIC
        ShpTemp
    );

    // Snh3
    MPtrs_[1].ref() = fSion
    (
        KaIN_,
        YPtrs_[10].internalField(), // SIN
        ShpTemp
    );

    // calc SohN
    // TODO: the original ADMno1 is quite inconsistent with the dimensions
    // TODO: maybe reverse the dimensionsScalar to scalar in para_?
    volScalarField::Internal SohN = KaW_ / ShpTemp;
    SohN.dimensions().reset(ShpTemp.dimensions()); 

    volScalarField::Internal E = 
    (
        IOPtrs_[0].internalField() - IOPtrs_[1].internalField() + ShpTemp 
      - SohN + (YPtrs_[10].internalField() - MPtrs_[1].internalField()) - EPtrs_[4] 
      - para_.MTOm() * (EPtrs_[3]/64.0 + EPtrs_[2]/112.0 + EPtrs_[1]/160.0 + EPtrs_[0]/208.0)
    );

    // DEBUG
    // Info<< ">>>\n" <<
    //        "E(x):\t" << max(E.field()) << "\n" << // endl;
    //        "x:\t" << max(ShpTemp.field()) << "\n" <<
    //        "SvaN:\t" << max(SvaN.field()) << "\n" <<
    //        "SbuN:\t" << max(SbuN.field()) << "\n" <<
    //        "SproN:\t" << max(SproN.field()) << "\n" <<
    //        "SacN:\t" << max(SacN.field()) << "\n" <<
    //        "Shco3N:\t" << max(Shco3N.field()) << "\n" <<
    //        "Snh3:\t" << max(MPtrs_[1].field()) << "\n" <<
    //        "SohN:\t" << max(SohN.field()) << "\n" << endl;

    return E;
}


volScalarField::Internal Foam::ADMno1::dfShp
(
    volScalarField::Internal& ShpTemp
)
{
    volScalarField::Internal dSvaN = dfSion
    (
        para_.Ka().va,
        YPtrs_[3].internalField(),
        ShpTemp
    );

    volScalarField::Internal dSbuN = dfSion
    (
        para_.Ka().bu,
        YPtrs_[4].internalField(),
        ShpTemp
    );

    volScalarField::Internal dSproN = dfSion
    (
        para_.Ka().pro,
        YPtrs_[5].internalField(),
        ShpTemp
    );

    volScalarField::Internal dSacN = dfSion
    (
        para_.Ka().ac,
        YPtrs_[6].internalField(),
        ShpTemp
    );

    volScalarField::Internal dShco3N = dfSion
    (
        Kaco2_,
        YPtrs_[9].internalField(), // SIC
        ShpTemp
    );

    volScalarField::Internal dSnh3 = dfSion // Snh3
    (
        KaIN_,
        YPtrs_[10].internalField(), // SIN
        ShpTemp
    );

    // calc SohN
    // TODO: the original ADMno1 is quite inconsistent with the dimensions
    // TODO: maybe reverse the dimensionsScalar to scalar in para_?
    volScalarField::Internal dSohN = - KaW_ / (ShpTemp * ShpTemp);
    dSohN.dimensions().reset(dSvaN.dimensions());

    dimensionedScalar uniField
    (
        dimless,
        One
    );

    return uniField - dSnh3 - dShco3N - dSohN  
         - para_.MTOm() * (dSacN/64.0 + dSproN/112.0 + dSbuN/160.0 + dSvaN/208.0);
}


void Foam::ADMno1::calcShp()
{
    //TODO: IO dictionary for these parameters
    scalar tol = 1e-12;
    label nIter = 1e3;
    label i = 0;

    // initial value of x, E and dEdx
    volScalarField::Internal x = ShP_;   // x = Shp
    volScalarField::Internal E = ShP_;   // E = dShp/dt
    volScalarField::Internal dE = ShP_;  // dE = (dShp/dt)/dShp
    
    do
    {
        E.field() = fShp(x).field();
        dE.field() = dfShp(x).field();
        x.field() = x.field() - E.field()/dE.field();
        i++;
    }
    while
    (
        // gMax(mag(E.field())) > tol &&
        max(mag(E.field())) > tol &&
        i < nIter
    );

    x.field() = min(max(x.field(), scalar(1e-16)), x.field());

    // if
    // (
    //     min(x.field()) < 0
    // )
    // {
    //     x.field() = 0.0*x.field() + 1e-16;
    // }

    Info<< "Newton-Raphson:\tSolving for Sh+" 
        << ", min Shp: " << min(x.field()) 
        << ", max Shp: " << max(x.field()) 
        << ", No Interations " << i << endl;

    // ShP
    ShP_ = x;
    pH_.field() = -log10(ShP_.field() / para_.MTOm());

    // update Sco2: Sco2 = SIC - Shco3N
    MPtrs_[0].ref() = YPtrs_[9] - EPtrs_[4];
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ADMno1::clear()
{

    forAll(dYPtrs_, i)
    {
        dYPtrs_[i] *= 0.0;
    }

    forAll(dGPtrs_, i)
    {
        dGPtrs_[i] *= 0.0;
    }
}

void Foam::ADMno1::correct
(
    const surfaceScalarField &flux,
    volScalarField& T
)
{
    //- Calculate thermal factor and adjust parameters
    calcThermal(T);

    //- Gas phase pressure
    gasPressure();

    //- Inhibition rates
    inhibitions();

    //- calculate raction rates
    kineticRate();

    //- calculate gas phase transfer rates
    gasPhaseRate();

    //- calculate dY with STOI
    dYUpdate(flux);

    //- calculate gas exit rates
    gasSourceRate();

    //- Acid-base calculations
    calcShp();

    //- Sh2 calculations
    calcSh2(flux);
}


tmp<fvScalarMatrix> Foam::ADMno1::R
(
    label i
) const
{
    DimensionedField<scalar, volMesh> dY = dYPtrs_[i];

        tmp<fvScalarMatrix> tSu
        (
            new fvScalarMatrix
            (
                YPtrs_[i],
                dY.dimensions()*dimVolume
                // dimMass/dimTime // <- for compressible flow and uses fvm::ddt(rho, Yi)
            )
        );

    fvScalarMatrix& Su = tSu.ref();
    
    // https:// /documentation/guides/latest/api/fvMatrix_8C_source.html#l01708
    Su += dY; 

    return tSu;
}; 


tmp<fvScalarMatrix> Foam::ADMno1::RG
(
    label i
) const
{
    DimensionedField<scalar, volMesh> dG = dGPtrs_[i];

        tmp<fvScalarMatrix> tSu
        (
            new fvScalarMatrix
            (
                GPtrs_[i],
                dG.dimensions()*dimVolume
                // dimMass/dimTime // <- for compressible flow and uses fvm::ddt(rho, Yi)
            )
        );

    fvScalarMatrix& Su = tSu.ref();
    
    // https:// /documentation/guides/latest/api/fvMatrix_8C_source.html#l01708
    Su += dG; 

    return tSu;
};

tmp<fvScalarMatrix> Foam::ADMno1::RIO
(
    label i
) const
{
    DimensionedField<scalar, volMesh> dIO = dIOPtrs_[i];

        tmp<fvScalarMatrix> tSu
        (
            new fvScalarMatrix
            (
                IOPtrs_[i],
                dIO.dimensions()*dimVolume
                // dimMass/dimTime // <- for compressible flow and uses fvm::ddt(rho, Yi)
            )
        );

    fvScalarMatrix& Su = tSu.ref();
    
    // https:// /documentation/guides/latest/api/fvMatrix_8C_source.html#l01708
    Su += dIO; 

    return tSu;
}; 

// ************************************************************************* //