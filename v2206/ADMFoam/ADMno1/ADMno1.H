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

    //- Total reaction tank volume [m3]
    const dimensionedScalar Vmesh_;

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

    //- Schmidt number
    scalar Sc_;

    //- Turbulent Schmidt number
    scalar Sct_;

    //- Dimensionless operating temperature field
    volScalarField TopDummy_;

    //- Thermal condition factor
    volScalarField::Internal fac_;

    //- Thermally adjusted parameters
    volScalarField::Internal KHh2_, KHch4_, KHco2_,
                             KaW_,  KaIN_,  Kaco2_;                    

    //- Ihibition factors
    dimensionedScalar nIaa_, nIac_, nIh2_;

    //- pH value
    volScalarField pH_;

    //- Concentration of proton ions (H+)
    volScalarField::Internal ShP_;

    //- List of pointers to the field of specie concentrations
    PtrList<volScalarField> YPtrs_;

    //- List of pointers to the field of specie concentrations 
    //  derivatives (source for local cell)
    PtrList<volScalarField::Internal> dYPtrs_;

    //- List of pointers to the field of gas concentrations
    PtrList<volScalarField> GPtrs_;

    //- List of pointers to the field of gas concentrations 
    //  derivatives (source for local cell)
    PtrList<volScalarField::Internal> dGPtrs_;

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

    //- List of pointers to gas trasfer rate
    PtrList<volScalarField::Internal> GRPtrs_;

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

    //- Cation concentration (inert)
    dimensionedScalar Scat_;

    //- Anion concentration (inert)
    dimensionedScalar San_;

    //- Time Scale
    dimensionedScalar tc_;

    //- Pirvate Helper Functions
    #include "ADMno1Helpers.H"
    

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

    //- Adjust for thermal condition
    void calcThermal(volScalarField& T);

    //- Gas phase pressure
    void gasPressure();

    //- Gas phase transfer rate
    void gasPhaseRate();

    //- Gas exit rates dG
    void gasSourceRate();

    //- Calculate inhibition rates
    void inhibitions();

    //- Calculate kinetic rates
    void kineticRate();

    //- Calculate derivative of species dY
    void dYUpdate(const surfaceScalarField &flux);

    //- Calculate acid-base reactions (iteratively)
    volScalarField::Internal fShp
    (
        volScalarField::Internal &ShpTemp
    );

    volScalarField::Internal dfShp
    (
        volScalarField::Internal &ShpTemp
    );

    void calcShp();

    //- Calculate Sh2 (iteratively)
    volScalarField::Internal fSh2
    (
        const surfaceScalarField &flux,
        volScalarField &Sh2Temp
    );

    volScalarField::Internal dfSh2
    (
        const surfaceScalarField &flux,
        volScalarField &Sh2Temp
    );

    void calcSh2
    (
        const surfaceScalarField &flux
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    ADMno1
    ( 
        const fvMesh& mesh,
        word runMode 
    );

    ADMno1
    ( 
        volScalarField& T,
        const fvMesh& mesh,
        const IOdictionary& ADMno1Dict
    );


    //- Selectors

    // template<class ADModel>
    // static autoPtr<ADModel> New (const fvMesh&);

    static autoPtr<ADMno1> New
    (
        volScalarField& T,
        const fvMesh&
    );


    //- Destructor

    ~ADMno1() = default;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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

    //- Return GPtrs_
    inline PtrList<volScalarField>& G();

    //- Return dGPtrs_
    inline PtrList<volScalarField::Internal>& dG();

    //- Return IOPtrs_
    inline PtrList<volScalarField>& IO();

    //- Return dIOPtrs_
    inline PtrList<volScalarField::Internal>& dIO();

    //- Clear derivative pointer lists
    void clear();

    //- Correct components consumption rate
    void correct(const surfaceScalarField &flux, volScalarField& T);

    //- Return components consumption rate matrix, i.e. source term.
    tmp<fvScalarMatrix> R(label i) const;

    //- Return components consumption rate matrix, i.e. source term.
    tmp<fvScalarMatrix> RG(label i) const;
    
    //- Return components consumption rate matrix, i.e. source term.
    tmp<fvScalarMatrix> RIO(label i) const;
   
    //- Return const access to the chemical source terms for specie, i
    // inline const volScalarField::Internal& RR ( const label i ) const;

    //- Return non const access to chemical source terms [kg/m3/s]
    // virtual volScalarField::Internal& RR ( const label i );

    // combustionModels/laminar.C -> template<class ReactionThermo>
                                  // Foam::tmp<Foam::fvScalarMatrix>
                                  // Foam::combustionModels::laminar<ReactionThermo>::R(volScalarField& Y) const


    void calcGasParameters()
    {
        Vfrac_.value() = this->get<scalar>("Vfrac");
        Vgas_ = Vmesh_ * Vfrac_ / (1 + Vfrac_);
        Vliq_ = Vmesh_ - Vgas_;

        // just a test too generate layers for case studies
        // Qin_.value() = 0.0;
        // Vgas_.value() = 300;
        // Vliq_.value() = 3400;
        // Vfrac_ = Vgas_/Vliq_;

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
            Info<< ">>> Benchmark case [Rose and Jeppsson, 2006]"<< endl;

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
        // KP_.value() = this->get<scalar>("Kpip") / 86400.0;
        
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