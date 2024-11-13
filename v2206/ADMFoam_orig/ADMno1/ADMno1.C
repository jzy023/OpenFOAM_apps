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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


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
        ADMno1Dict.lookupOrDefault("qin", 0.00)
        // ADMno1Dict.lookupOrDefault("qin", 178.4674) // benchmark 
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
    para_(ADMno1Dict.get<word>("mode")),
    Sc_(ADMno1Dict.lookupOrDefault("Sc", 1.0)),
    R_(ADMno1Dict.lookupOrDefault("R", 0.083145)),
    KP_(ADMno1Dict.lookupOrDefault("Kpip", 5e4)),
    // Vfrac_(ADMno1Dict.lookupOrDefault("Vfrac", 0.0294118)), // 100/3400
    Vfrac_(ADMno1Dict.lookupOrDefault("Vfrac", 0.0882353)), // 300/3400
    Pext_
    (
        "Pext", 
        dimPressure,
        ADMno1Dict.lookupOrDefault("Pext", 1.013)
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
        ADMno1Dict.lookupOrDefault("San", 0.0052)
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

    for(label i = 0; i < nSpecies; i++)
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

    for(label i = 0; i < namesGaseous.size(); i++)
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
                ADMno1Dict.lookupOrDefault("San", 0.0052)
            )
        )
    );

    
    //- Initializing derivatives

    dIOPtrs_.resize(namesIons.size());

    for(label i = 0; i < namesIons.size(); i++)
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
    MPtrs_[0].ref() = YPtrs_[9] - EPtrs_[4]; // Sco2 = SIC - Shco3N

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

    fac_ = (1.0 / para_.Tbase().value() - 1.0 / TopDummy_) / (100.0 * R_);
    
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


void Foam::ADMno1::calcTC()
{

    forAll(dYPtrs_, i)
    {
        if (i == 7)
        {
            continue;
        }
        
        tc_ = min
        (
            min(mag(YPtrs_[i]/dYPtrs_[i])), 
            tc_
        );
    }
    
    Info<< "time scale: " << tc_.value() << endl;
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

    // //- 
    // // calcTC();

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