    /*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |  
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C)  Jeremy Z. Yan
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

Class
    ADMno1
    >>> inspired from reactingFoam/CombustionModel<psiReactionThermo>

Description
    Anaerobic Digestion ADModel No.1 class.


\*---------------------------------------------------------------------------*/

#ifndef ADMno1_H
#define ADMno1_H

#include "ADMno1Parameter.H"
#include "autoPtr.H"
#include "PtrListDictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class ADMno1
:
    public IOdictionary
{

private:
    //- Reference of VoF of liquid phase from mixture 
    const volScalarField& alpha1_;

    //- Total reaction tank volume [m3]
    const dimensionedScalar Vmesh_;

    // ------------------------- Gas Phase Parameters ---------------------- //
    //- Volume content [m3 s-1]
    //  benchmark purpose only
    dimensionedScalar Qin_; 

    //- Total gas phase volume [m3]
    dimensionedScalar Vgas_;

    //- Total liquid phase volume [m3]
    dimensionedScalar Vliq_;

    //- Gas-Liquid volume fraction in tank 
    dimensionedScalar Vfrac_;
    
    //- Dimensionless gas outflow rate per pascal of presure difference
    dimensionedScalar qGas_;

    //- Pipe resistance coefficient [m³s⁻¹Bar⁻¹]
    dimensionedScalar KP_;

    // TODO: [s-1]?
    // Gas transfer rate coefficients [d-1]
    dimensionedScalar kLa_;

    //- Vapour pressure of water [Pa]
    dimensionedScalar Pvap_; 

    //- Gas exit ambient pressure
    dimensionedScalar Pext_;

    //- Gas exit ambient pressure
    volScalarField Pgas_; 

    //- Gas constant
    // dimensionedScalar R_;
    scalar R_;

    // -------------------------- General Parameters ----------------------- //
    //- Operating under benchmark mode or not
    const bool isBenchmark_;

    //- ADM no1 operation mode
    const word runMode_;

    //- ADM no1 parameters
    admPara para_;

    //- Schmidt and turbulent Schmidt number
    scalar Sc_, Sct_;

    // -------------------------- Thermal Parameters ----------------------- //
    // TODO: fix dimension and move to [const volScalarField& Top]
    //- Dimensionless operating temperature field
    volScalarField TopDummy_;

    //- Average Temperature for gaseous partial pressure
    dimensionedScalar TopAve_;

    //- Thermal condition factor
    volScalarField::Internal fac_;

    //- Thermally adjusted parameters
    volScalarField::Internal KaW_,  KaIN_,  Kaco2_;                    
    dimensionedScalar        KHh2_, KHch4_, KHco2_;

    //- Ihibition factors
    dimensionedScalar nIaa_, nIac_, nIh2_;

    // ------------------------ Species Concentrations --------------------- //
    //- Concentration of proton ions (H+)
    volScalarField::Internal ShP_;

    //- pH value
    volScalarField pH_;

    //- List of pointers to the field of specie concentrations
    PtrList<volScalarField> YPtrs_;

    //- List of pointers to the field of specie concentrations 
    //  derivatives (source for local cell)
    PtrList<volScalarField::Internal> dYPtrs_;

    // Averaged concentrations of gaseous species (perfect mixture)
    PtrList<dimensionedScalar> GAvePtrs_;

    // Averaged volumetric mass transfer rate of gaseous species (perfect mixture)
    PtrList<dimensionedScalar> GRAvePtrs_;

    //- List of pointers to the field of gas concentrations derivatives
    PtrList<dimensionedScalar> dGAvePtrs_;

    //- List of pointers to the field of medians concentrations
    PtrList<volScalarField> MPtrs_;

    //- List of pointers to the field of ions concentrations
    PtrList<volScalarField> IOPtrs_;

    //- List of pointers to the field of ion concentrations derivatives
    PtrList<volScalarField::Internal> dIOPtrs_;

    //- List of pointers to the field of electrolytes concentrations
    PtrList<volScalarField::Internal> EPtrs_;

    //- List of pointers to inhibition 
    PtrList<volScalarField::Internal> IPtrs_;

    //- List of pointers to kinetic rate per specie [kg/m3/s] and [mol/m3/s]
    PtrList<volScalarField::Internal> KRPtrs_;
    
    //- Cation and Anion concentration (inert)
    dimensionedScalar Scat_, San_;

    // --------------------------------------------------------------------- //
    //- Pirvate Helper Functions
    #include "ADMno1Helpers.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

    //- Adjust for thermal condition
    void calcThermal(const volScalarField& T);

    //- Gas phase pressure
    void gasPressure();

    //- Gas phase transfer rate
    // void gasPhaseRate();
    void gasPhaseRate(const volScalarField& kLaCells);

    //- Gas exit rates dG
    void gasSourceRate();
    // void gasSourceRate(const volScalarField& actCells);

    //- Calculate inhibition rates
    void inhibitions();

    //- Calculate kinetic rates
    void kineticRate();

    //- Calculate derivative of species dY
    void dYUpdate(const surfaceScalarField& phi);

    //- Calculate acid-base reactions (Newton-Raphson)
    volScalarField::Internal fShp
    (
        volScalarField::Internal &ShpTemp
    );

    volScalarField::Internal dfShp
    (
        volScalarField::Internal &ShpTemp
    );

    void calcShp();

    //- Calculate Sh2 (Newton-Raphson)
    volScalarField::Internal fSh2
    (
        const volScalarField& kLaCells,
        const surfaceScalarField& phi,
        volScalarField& Sh2Temp
    );

    volScalarField::Internal dfSh2
    (
        const volScalarField& kLaCells,
        const surfaceScalarField& phi,
        volScalarField& Sh2Temp
    );

    void calcSh2
    (
        const volScalarField& kLaCells,
        const surfaceScalarField& phi
    );
    

public:
 
    // //- Run time type information
    // //> runTimeSelection mechanism: 
    // //> https://openfoamwiki.net/index.php/OpenFOAM_guide/runTimeSelection_mechanism
    // TypeName("ADMno1");

    // // Declare run-time constructor selection table
    // declareRunTimeSelectionTable
    // (
    //     autoPtr,
    //     ADMno1,
    //     fvMesh,
    //     (
    //         volScalarField& T,
    //         const fvMesh& mesh,
    //         word runMode
    //     ),
    //     (
    //         T,
    //         mesh,
    //         runMode
    //     )
    // );

    //- Default name of the ADMno1 properties dictionary
    static const word propertiesName;

    //- Species names for soluable
    //- !!! Sh2 is calculated iteratively 
    const wordList namesSoluable = 
    //      0      1      2      3      4      5
          {"Ssu", "Saa", "Sfa", "Sva", "Sbu", "Spro",
           "Sac", "Sh2", "Sch4", "SIC", "SIN", "SI"};
    //      6      7      8       9      10     11
  
    //- Species names for particulate
    const wordList namesParticulate = 
    //    12    13     14     15     16     17
          {"Xc", "Xch", "Xpr", "Xli", "Xsu", "Xaa",
           "Xfa", "Xc4", "Xpro", "Xac", "Xh2", "XI"};
    //      18     19     20      21     22     23

    //- Species names for gaseouse
    const wordList namesGaseous = 
          {"Gh2", "Gch4", "Gco2"};

    //- Species names for electrolytes
    const wordList namesElectrolytes = 
          {"SvaN", "SbuN", "SproN", "SacN", "Shco3N"};
    //      0       1       2        3       4

    //- Calculated by simple subtractions
    const wordList namesMedians = 
          {"Sco2", "Snh3"};

    //- Ions substances
    const wordList namesIons = 
          {"Scat", "San"};

    //- Default gas concentration [Rose and Jeppsson, 2006]
    const scalarList benchmarkGaseous = 
    {1.1032e-5, 1.6535, 13.5};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    ADMno1
    ( 
        const fvMesh& mesh,
        word runMode 
    );

    ADMno1
    ( 
        const volScalarField& alpha1,
        const volScalarField& T,
        const fvMesh& mesh,
        const IOdictionary& ADMno1Dict
    );


    //- Selectors

    // template<class ADModel>
    // static autoPtr<ADModel> New (const fvMesh&);

    static autoPtr<ADMno1> New
    (
        const volScalarField& alpha1,
        const volScalarField& T,
        const fvMesh&
    );


    //- Destructor

    ~ADMno1() = default;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
    
    //- Initialization
    void init(const volScalarField& T);
    
    //- Clear derivative pointer lists
    void clear();

    //- Correct components consumption rate
    void correct
    (
        const volScalarField& kLaCells,
        const surfaceScalarField& phi,
        const volScalarField& T
    );

    //- Inline Member Functions
    //- Return Schmidt number
    inline scalar Sc();

    //- Return turbulent Schmidt number
    inline scalar Sct();

    //- Return a pointer of admno1 parameters
    inline admPara* paraPtr();
    
    //- Return YPtrs_
    inline PtrList<volScalarField>& Y();

    //- Return dYPtrs_
    inline PtrList<volScalarField::Internal>& dY();

    //- Return GAvePtrs_
    inline PtrList<dimensionedScalar>& GAve();

    //- Return GRAvePtrs_
    inline PtrList<dimensionedScalar>& GRAve();

    //- Return dGAvePtrs_
    inline PtrList<dimensionedScalar>& dGAve();

    //- Return MPtrs_
    inline PtrList<volScalarField>& M();

    //- Return IOPtrs_
    inline PtrList<volScalarField>& IO();

    //- Return dIOPtrs_
    inline PtrList<volScalarField::Internal>& dIO();

    // //- Return components consumption rate matrix, i.e. source term.
    // tmp<fvScalarMatrix> R(label i) const;

    // //- Return components consumption rate matrix, i.e. source term.
    // tmp<fvScalarMatrix> RG(label i) const;
    
    // //- Return components consumption rate matrix, i.e. source term.
    // tmp<fvScalarMatrix> RIO(label i) const;
   
    //- Return const access to the chemical source terms for specie, i
    // inline const volScalarField::Internal& RR ( const label i ) const;

    //- Return non const access to chemical source terms [kg/m3/s]
    // virtual volScalarField::Internal& RR ( const label i );

    // combustionModels/laminar.C -> template<class ReactionThermo>
                                  // Foam::tmp<Foam::fvScalarMatrix>
                                  // Foam::combustionModels::laminar<ReactionThermo>::R(volScalarField& Y) const
                                  

    // testing ------------------------------------------------------------
    //- Return Henry's constant
    // TODO: maybe return PtrList<scalar> instead?
    // inline volScalarField::Internal& KHh2()
    // {
    //     return KHh2_;
    // };

    // inline volScalarField::Internal& KHch4()
    // {
    //     return KHch4_;
    // };

    // inline volScalarField::Internal& KHco2()
    // {
    //     return KHco2_;
    // };

    // inline const dimensionedScalar PgasAve()
    // {
    //     return Pgas_.weightedAverage(Pgas_.mesh().V());
    // };

    void calcGasParameters()
    {
        // calculate with alpha1
        dimensionedScalar alphaLiq = alpha1_.weightedAverage(alpha1_.mesh().V());
        Vfrac_ = 1.0/alphaLiq - 1.0;
        Vliq_ = alphaLiq * Vmesh_;
        Vgas_ = Vmesh_ - Vliq_;
        qGas_.value() = KP_.value() / Vgas_.value();
    }

    word checkBenchmark()
    {
        // initialize dimensions
        Qin_.dimensions().reset(dimVolume/dimTime);
        Vgas_.dimensions().reset(dimVolume);
        Vliq_.dimensions().reset(dimVolume);
        Vfrac_.dimensions().reset(dimless);
        qGas_.dimensions().reset(dimless/dimPressure/dimTime);
        KP_.dimensions().reset(dimVolume/dimPressure/dimTime);

        // >>> if benchmark case
        if (isBenchmark_)
        {
            Info<< ">>> Benchmark case [Rose and Jeppsson, 2006]\n"<< endl;

            Qin_.value() = 178.4674;
            Vgas_.value() = 300;
            Vliq_.value() = 3400;
            Vfrac_ = Vgas_/Vliq_;

            // KP: [m3 * s-1 * Pa-1]
            // KP_.value() = para_.DTOS() * 5e4 / para_.BTOP();
            KP_.value() = 5e4 / (1e5 * 86400.0);
            qGas_ = KP_ / Vgas_;

            return "MesoSolid";
        }

        // >>> if not benchmark case
        Qin_.value() = 0.0;
        KP_.value() = this->get<scalar>("Kpip");

        // calculate with alpha1
        calcGasParameters();

        return this->get<word>("mode");;
    }

    
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#include "ADMno1I.H"
#ifdef NoRepository
    #include "ADMno1Templates.C"
#endif
 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //