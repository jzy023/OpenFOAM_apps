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
    const volScalarField& alpha1,
    const volScalarField& T,
    const fvMesh& mesh,
    const IOdictionary& ADMno1Dict
)
:
    IOdictionary(ADMno1Dict),
    alpha1_
    (
        alpha1
    ),
    Vmesh_
    (
        dimVolume,
        gSum(alpha1_.mesh().V())
    ),
    kLa_
    (
        "kLa",
        dimless/dimTime,
        this->lookupOrDefault
        (
            "kLa",
            200
        ) 
    ),
    kLaCellsAve_
    (
        "kLaCellsAve",
        dimless/dimTime,
        Zero
    ),
    Pvap_
    (
        "Pvap", 
        dimPressure,
        this->get<scalar>("Pvap")
    ),
    Pext_
    (
        "Pext", 
        dimPressure,
        this->lookupOrDefault
        ("Pext", 1e5 * 1.013)
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
    R_
    (
        this->lookupOrDefault("R", 8.3145)
    ),
    // =============================================================
    isBenchmark_
    (
        this->lookupOrDefault("benchmark", false)
    ),
    runMode_
    (
        checkBenchmark()
    ),
    para_
    (
        runMode_
    ),
    Sc_
    (
        this->lookupOrDefault("Sc", 1.0)
    ),
    Sct_
    (
        this->lookupOrDefault("Sct", 0.2)
    ),
    TopDummy_(T),
    TopAve_(308.15),
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
            this->lookupOrDefault("fac", 1)
        )
    ),
    // KHh2_(fac_),
    // KHch4_(fac_),
    // KHco2_(fac_),
    KaW_(fac_),
    KaIN_(fac_),                    
    Kaco2_(fac_),
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
            this->lookupOrDefault("pH", 7.26)
        )
    ),
    Scat_
    (
        "Scat",
        dimMass/dimVolume,
        this->lookupOrDefault("Scat", 0.00)
    ),
    San_
    (
        "San",
        dimMass/dimVolume,
        this->lookupOrDefault("San", 0.0052 * para_.MTOm())
    )
{

    Info<< "\nSelecting ADM no1 operation mode " << this->get<word>("mode") << endl;

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
                    IOobject::groupName(namesSoluable[i], "ADM"),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    // IOobject::NO_WRITE
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    dimMass/dimVolume, 
                    Zero
                )
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
                    IOobject::groupName(namesParticulate[i], "ADM"),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    dimMass/dimVolume, 
                    Zero
                )
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

    GAvePtrs_.resize(namesGaseous.size());

    forAll(namesGaseous, i)
    {
        GAvePtrs_.set
        (
            i,
            new dimensionedScalar
            (
                namesGaseous[i], 
                YPtrs_[0].dimensions(),
                this->lookupOrDefault
                (
                    namesGaseous[i],
                    benchmarkGaseous[i] // TODO: move to checkBenchmark()
                )
            )
        );
    }
    
    //- Initializing derivatives

    dGAvePtrs_.resize(namesGaseous.size());

    for (label i = 0; i < namesGaseous.size(); i++)
    {
        dGAvePtrs_.set
        (
            i,
            new dimensionedScalar
            (
                "d" + namesGaseous[i],
                GAvePtrs_[0].dimensions()/dimTime, 
                Zero
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
                    IOobject::groupName(namesMedians[i], "ADM"),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
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
                this->lookupOrDefault("Scat", 0.00)
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
                this->lookupOrDefault("San", 0.0052)
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
                    IOobject::groupName("Inh", Foam::name(i)),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
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

    GRPtrs_.resize(namesGaseous.size());

    forAll(namesGaseous, i)
    {
        GRPtrs_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    namesGaseous[i] + ".R",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    GAvePtrs_[0].dimensions()/dimTime,
                    Zero
                )
            )
        );
    }
    
    GRAvePtrs_.resize(namesGaseous.size());

    forAll(namesGaseous, i)
    {
        GRAvePtrs_.set
        (
            i,
            new dimensionedScalar
            (
                namesGaseous[i] + "Ave.R",
                GAvePtrs_[0].dimensions()/dimTime,
                Zero
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // reset dimensions 
    para_.setParaDim(YPtrs_[0].dimensions());

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

    // YPtrs_[7].field() = 2.5055e-7;

    // DEBUG
    Info<<   ">>> is benchmark case: " << isBenchmark_
        << "\n>>> Qin_: "   << Qin_.value()
        << "\n>>> Vgas_: "  << Vgas_.value()
        << "\n>>> Vliq_: "  << Vliq_.value()
        << "\n>>> Vfrac_: " << Vfrac_.value()
        << "\n>>> qGas_: "  << qGas_.value() << endl;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
 
Foam::autoPtr<Foam::ADMno1> Foam::ADMno1::New
(
    const volScalarField& alpha1,
    const volScalarField& T,
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
    ADMno1* reactionPtr = new ADMno1
    (
        alpha1,
        T,
        mesh,
        ADMno1Dict
    );
    return autoPtr<ADMno1>(reactionPtr);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ADMno1::calcThermal  
(
    const volScalarField& T
)
{
    // DEBUG MULTI
    TopDummy_.field() = T.field();
    TopAve_ = TopDummy_.weightedAverage(YPtrs_[0].mesh().V());

    fac_ = (1.0 / para_.Tbase().value() - 1.0 / TopDummy_) / R_;
    
    Kaco2_ = para_.Ka().co2 * exp(7646.0 * fac_);
    KaIN_ = para_.Ka().IN * exp(51965.0 * fac_);
    KaW_ = para_.Ka().W * exp(55900.0 * fac_);

    // Average of field for gas transport
    dimensionedScalar facAve =
    (
        (1.0 / para_.Tbase().value() - 1.0 / TopAve_.value()) / R_
    );
    
    KHh2_ = para_.KH().h2 * exp(-4180.0 * facAve);
    KHch4_ = para_.KH().ch4 * exp(-14240.0 * facAve);
    KHco2_ = para_.KH().co2 * exp(-19410.0 * facAve);
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
            YPtrs_[0].dimensions(),
            Zero
        )
    );

    Ph2o.field() = Pvap_ * exp(5290.0 * fac_ * R_);
    
    Pgas_.field() = 
    (
        Ph2o + R_ * TopDummy_
     * (para_.MTOm() * GAvePtrs_[0] / 16.0 + para_.MTOm() * GAvePtrs_[1] / 64.0 + GAvePtrs_[2])
    );
}


void Foam::ADMno1::gasPhaseRate
(
    const volScalarField& kLaCells
)
{
    // New method 1 ---------------------------------------------------------------------------
    dimensionedScalar Sh2Ave = YPtrs_[7].weightedAverage(YPtrs_[0].mesh().V());
    dimensionedScalar Sch4Ave = YPtrs_[8].weightedAverage(YPtrs_[0].mesh().V());
    dimensionedScalar Sco2Ave = MPtrs_[0].weightedAverage(YPtrs_[0].mesh().V());

    Info<< ">>> Sh2Ave = " << Sh2Ave.value() << endl;
    Info<< ">>> Sch4Ave = " << Sch4Ave.value() << endl;
    Info<< ">>> Sco2Ave = " << Sco2Ave.value() << endl;
    
    GRPtrs_[0] = // <-- kg COD m-3 s-1
    (   // kLa_new = [some correction factor] * isCellInterface_ + alphaW_ * isCellsWall_
        kLaCells * (Sh2Ave - R_ * TopAve_ * GAvePtrs_[0] * KHh2_)
    );

    GRPtrs_[1] = // <-- kg COD m-3 s-1
    (   // kLa_new = [some correction factor] * isCellInterface_ + alphaW_ * isCellsWall_
        kLaCells * (Sch4Ave - R_ * TopAve_ * GAvePtrs_[1] * KHch4_)
    );

    GRPtrs_[2] = // <-- mol COD m-3 s-1
    (   // Sco2 instead of SIC
        // kLa_new = [some correction factor] * isCellInterface_ + alphaW_ * isCellsWall_ 
        kLaCells * (Sco2Ave - R_ * TopAve_ * GAvePtrs_[2] * KHco2_)
    );

    GRAvePtrs_[0] = GRPtrs_[0].weightedAverage(GRPtrs_[0].mesh().V());
    GRAvePtrs_[1] = GRPtrs_[1].weightedAverage(GRPtrs_[1].mesh().V());
    GRAvePtrs_[2] = GRPtrs_[2].weightedAverage(GRPtrs_[2].mesh().V());

    // kLaCellsAve_ = kLaCells.weightedAverage(kLaCells.mesh().V());
    // kLaCellsAve_ = para_.DTOS() * kLa_;

    // GRAvePtrs_[0] = // <-- kg COD m-3 s-1
    // (   // kLa_new = [some correction factor] * isCellInterface_ + alphaW_ * isCellsWall_ 
    //     kLaCellsAve_ * (Sh2Ave - R_ * TopAve_ * GAvePtrs_[0] * KHh2_)
    // );

    // GRAvePtrs_[1] = // <-- kg COD m-3 s-1
    // (   // kLa_new = [some correction factor] * isCellInterface_ + alphaW_ * isCellsWall_ 
    //     kLaCellsAve_ * (Sch4Ave - R_ * TopAve_ * GAvePtrs_[1] * KHch4_)
    // );

    // GRAvePtrs_[2] = // <-- mol COD m-3 s-1
    // (   // Sco2 instead of SIC
    //     // kLa_new = [some correction factor] * isCellInterface_ + alphaW_ * isCellsWall_  
    //     kLaCellsAve_ * (Sco2Ave - R_ * TopAve_ * GAvePtrs_[2] * KHco2_)
    // );

    // DEBUG
    // Info<< ">>> kLa [s^-1] ADMno1: " << (para_.DTOS() * kLa_).value() << "\n"
    //     << ">>> kLa [s^-1] ADMno1Multi: " << kLaCellsAve_.value() << endl;

    forAll(GAvePtrs_, i)
    {
        Info<< ">>> "                 << GAvePtrs_[i].name()
            << " concentration = "    << GAvePtrs_[i].value()
            << ", generation rate = " << GRAvePtrs_[i].value() << endl;
    }
}

// Benchmark
// force kLaCells to a fixed value regardless of the multiphase model
void Foam::ADMno1::gasPhaseRateBenchmark()
{
    const dimensionedScalar kLaCells = para_.DTOS() * kLa_;

    dimensionedScalar Sh2Ave = YPtrs_[7].weightedAverage(YPtrs_[0].mesh().V());
    dimensionedScalar Sch4Ave = YPtrs_[8].weightedAverage(YPtrs_[0].mesh().V());
    dimensionedScalar Sco2Ave = MPtrs_[0].weightedAverage(YPtrs_[0].mesh().V());

    Info<< ">>> Sh2Ave = " << Sh2Ave.value() << endl;
    Info<< ">>> Sch4Ave = " << Sch4Ave.value() << endl;
    Info<< ">>> Sco2Ave = " << Sco2Ave.value() << endl;
    
    GRPtrs_[0] = // <-- kg COD m-3 s-1
    (   // kLa_new = [some correction factor] * isCellInterface_ + alphaW_ * isCellsWall_
        kLaCells * (Sh2Ave - R_ * TopAve_ * GAvePtrs_[0] * KHh2_)
    );

    GRPtrs_[1] = // <-- kg COD m-3 s-1
    (   // kLa_new = [some correction factor] * isCellInterface_ + alphaW_ * isCellsWall_
        kLaCells * (Sch4Ave - R_ * TopAve_ * GAvePtrs_[1] * KHch4_)
    );

    GRPtrs_[2] = // <-- mol COD m-3 s-1
    (   // Sco2 instead of SIC
        // kLa_new = [some correction factor] * isCellInterface_ + alphaW_ * isCellsWall_ 
        kLaCells * (Sco2Ave - R_ * TopAve_ * GAvePtrs_[2] * KHco2_)
    );

    GRAvePtrs_[0] = GRPtrs_[0].weightedAverage(GRPtrs_[0].mesh().V());
    GRAvePtrs_[1] = GRPtrs_[1].weightedAverage(GRPtrs_[1].mesh().V());
    GRAvePtrs_[2] = GRPtrs_[2].weightedAverage(GRPtrs_[2].mesh().V());

    forAll(GAvePtrs_, i)
    {
        Info<< ">>> "                 << GAvePtrs_[i].name()
            << " concentration = "    << GAvePtrs_[i].value()
            << ", generation rate = " << GRAvePtrs_[i].value() << endl;
    }
}


void Foam::ADMno1::gasSourceRate()
{
    volScalarField qGasField = qGas_ * (Pgas_ - Pext_);
    dimensionedScalar qGasAve = qGasField.weightedAverage(qGasField.mesh().V());

    // New method 1 ---------------------------------------------------------------------------
    dGAvePtrs_[0] = 
    (
        (GRAvePtrs_[0] * (1. / Vfrac_))    // <-- multiphase mass transfer
      - (GAvePtrs_[0] * qGasAve) // <-- gas released from chamber
    );

    dGAvePtrs_[1] = 
    (
        (GRAvePtrs_[1] * (1. / Vfrac_))    // <-- multiphase mass transfer
      - (GAvePtrs_[1] * qGasAve) // <-- gas released from chamber
    );

    dGAvePtrs_[2] = 
    (
        (GRAvePtrs_[2] * (1. / Vfrac_))    // <-- multiphase mass transfer
      - (GAvePtrs_[2] * qGasAve) // <-- gas released from chamber
    );
};


//- Functions for Sh2 calculations
volScalarField::Internal Foam::ADMno1::fSh2
(
    const surfaceScalarField& phi,
    volScalarField& Sh2Temp
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

    volScalarField conv
    (
        IOobject
        (
            "convSh2",
            alpha1_.mesh().time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar
        (
            dimDensity/dimTime, 
            Zero
        )
    );

    if (isBenchmark_)
    {
        conv = 
        (
            para_.DTOS() * (Qin_/Vliq_) * (para_.INFLOW(7) - Sh2Temp)
        //   + fvc::div(phi, Sh2Temp, "div(phi,Yi)") // ?
        );
    }
    else
    {
        // TODO: this div needs to be changed from admMixture->div() 
        //       to account for interface-corrected mass flux
        // conv = fvc::div(phi, Sh2Temp);
    }

    dimensionedScalar Sh2Ave = Sh2Temp.weightedAverage(YPtrs_[0].mesh().V());

    dimensionedScalar GRSh2Temp = 
    (
        kLaCellsAve_
      * (Sh2Ave - R_ * TopAve_ * GAvePtrs_[0] * KHh2_)
    );

    //     reaction + convection - fGasRhoH2(paraPtr, Sh2);
    return concPerComponent(7, KRPtrs_temp) + conv - GRSh2Temp;
}


volScalarField::Internal Foam::ADMno1::dfSh2
(
    const surfaceScalarField& phi,
    volScalarField& Sh2Temp
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

    volScalarField dConv
    (
        IOobject
        (
            "dConvSh2",
            alpha1_.mesh().time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar
        (
            dimless/dimTime, 
            Zero
        )
    );

    if (isBenchmark_)
    {
        dConv = 
        (
           - para_.DTOS() * (Qin_/Vliq_)
        //    + fvc::div(phi)
        );
    }
    else
    {
        // TODO: this div needs to be changed from admMixture->div() 
        //       to account for interface-corrected mass flux
        // dConv = fvc::div(phi);
    }
    
    dimensionedScalar dGRSh2Temp = para_.DTOS() * kLa_;
    // volScalarField& dGRSh2Temp = kLaCellsAve_;
    // dimensionedScalar dGRSh2Temp = kLaCellsAve_;

    //     dReaction + dConvection - dfGasRhoH2(paraPtr, Sh2);
    return concPerComponent(7, dKRPtrs_temp) + dConv - dGRSh2Temp;
}


void Foam::ADMno1::calcSh2
(
    const surfaceScalarField& phi
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
        E.field() = fSh2(phi, x).field();
        dE.field() = dfSh2(phi, x).field();
        x.field() = x.field() - E.field()/dE.field();
        i++;
    }
    while
    (
        gMax(mag(E.field())) > tol &&
        i < nIter
    );

    // safe guard 
    scalar range = 0.1;
    dimensionedScalar xAve = x.weightedAverage(x.mesh().V());
    x.field() = min
    (
        max
        (
            x.field(), xAve.value() * (1 - range)
        ), 
        xAve.value() * (1 + range)
    );

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
    const surfaceScalarField& phi
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
    dYPtrs_[8] -= GRAvePtrs_[1]; // Sch4 - Gch4
    dYPtrs_[9] -= GRAvePtrs_[2]; // SIC - Gco2
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
    volScalarField::Internal SohN = KaW_ / ShpTemp;
    SohN.dimensions().reset(ShpTemp.dimensions()); 

    volScalarField::Internal E = 
    (
        IOPtrs_[0].internalField() - IOPtrs_[1].internalField() + ShpTemp 
      - SohN + (YPtrs_[10].internalField() - MPtrs_[1].internalField()) - EPtrs_[4] 
      - para_.MTOm() * (EPtrs_[3]/64.0 + EPtrs_[2]/112.0 + EPtrs_[1]/160.0 + EPtrs_[0]/208.0)
    );

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
        gMax(mag(E.field())) > tol &&
        i < nIter
    );

    // safe guard 
    scalar range = 2.0;
    dimensionedScalar xAve = x.weightedAverage(x.mesh().V());
    x.field() = min
    (
        max
        (
            x.field(), xAve.value() / range
        ), 
        xAve.value() * range
    );

    // DEBUG
    Info<< ">>> ShpAve = " << x.weightedAverage(x.mesh().V()).value() << endl;
    //

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

void Foam::ADMno1::init
(
    const volScalarField& T
)
{
    nIaa_ = 3.0 / (para_.pHL().ULaa - para_.pHL().LLaa);  // aa
    nIac_ = 3.0 / (para_.pHL().ULac - para_.pHL().LLac);  // ac
    nIh2_ = 3.0 / (para_.pHL().ULh2 - para_.pHL().LLh2);  // h2
    
    calcThermal(T);
    
    MPtrs_[0].ref() = YPtrs_[9] - EPtrs_[4]; // Sco2 = SIC - Shco3N
}


void Foam::ADMno1::clear()
{

    forAll(dYPtrs_, i)
    {
        dYPtrs_[i] *= 0.0;
    }

    forAll(dGAvePtrs_, i)
    {
        dGAvePtrs_[i] *= 0.0;
    }
}


void Foam::ADMno1::correct
(
    const volScalarField& kLaCells,
    const surfaceScalarField& phi,
    const volScalarField& T
)
{
    //- Calculate gas transport parameters
    if (!isBenchmark_)
    {
        calcGasParameters();
    }
    
    //- Calculate thermal factor and adjust parameters
    calcThermal(T);

    //- Gas phase pressure
    gasPressure();

    //- Inhibition rates
    inhibitions();

    //- calculate raction rates
    kineticRate();

    //- calculate gas phase transfer rates
    gasPhaseRate(kLaCells);
    // gasPhaseRateBenchmark();

    //- calculate dY with STOI
    dYUpdate(phi);

    //- calculate gas exit rates
    gasSourceRate();

    //- Acid-base calculations
    calcShp();

    //- Sh2 calculations
    calcSh2(phi);
}


// tmp<fvScalarMatrix> Foam::ADMno1::R
// (
//     label i
// ) const
// {
//     DimensionedField<scalar, volMesh> dY = dYPtrs_[i];

//         tmp<fvScalarMatrix> tSu
//         (
//             new fvScalarMatrix
//             (
//                 YPtrs_[i],
//                 dY.dimensions()*dimVolume
//                 // dimMass/dimTime // <- for compressible flow and uses fvm::ddt(rho, Yi)
//             )
//         );

//     fvScalarMatrix& Su = tSu.ref();
    
//     // https:// /documentation/guides/latest/api/fvMatrix_8C_source.html#l01708
//     Su += dY; 

//     return tSu;
// }; 


// tmp<fvScalarMatrix> Foam::ADMno1::RG
// (
//     label i
// ) const
// {
//     DimensionedField<scalar, volMesh> dG = dGPtrs_[i];

//         tmp<fvScalarMatrix> tSu
//         (
//             new fvScalarMatrix
//             (
//                 GPtrs_[i],
//                 dG.dimensions()*dimVolume
//                 // dimMass/dimTime // <- for compressible flow and uses fvm::ddt(rho, Yi)
//             )
//         );

//     fvScalarMatrix& Su = tSu.ref();
    
//     // https:// /documentation/guides/latest/api/fvMatrix_8C_source.html#l01708
//     Su += dG; 

//     return tSu;
// };

// tmp<fvScalarMatrix> Foam::ADMno1::RIO
// (
//     label i
// ) const
// {
//     DimensionedField<scalar, volMesh> dIO = dIOPtrs_[i];

//         tmp<fvScalarMatrix> tSu
//         (
//             new fvScalarMatrix
//             (
//                 IOPtrs_[i],
//                 dIO.dimensions()*dimVolume
//                 // dimMass/dimTime // <- for compressible flow and uses fvm::ddt(rho, Yi)
//             )
//         );

//     fvScalarMatrix& Su = tSu.ref();
    
//     // https:// /documentation/guides/latest/api/fvMatrix_8C_source.html#l01708
//     Su += dIO; 

//     return tSu;
// }; 

// ************************************************************************* //