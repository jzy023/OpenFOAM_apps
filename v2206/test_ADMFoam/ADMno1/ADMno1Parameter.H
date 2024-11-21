/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |  
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    ADM::ADMno1Parameter

Description
    H file for declairing AMDno1 parameter class

\*---------------------------------------------------------------------------*/

#ifndef ADMno1Parameter_H
#define ADMno1Parameter_H

#include "ADMno1Types.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class admPara
{
private:

    //- Operation modes (thermal conditions for biomass)
    const wordList namesOpMode = {"Meso", "MesoSolid", "Thermo"}; 

    //- Error message
    int errMessage;

    //- Day-to-second conversion
    const scalar DTOS_ = 1.0/86400.0;

    //- Base dimension (kg/m3; mol/m3; mass fraction)
    dimensionSet ds_;

    //- Base operational temperature
    dimensionedScalar Tbase_;   

    //- Operational temperature
    dimensionedScalar Top_;   

    //- Decay rate [1/day]
    decayRate kDec_; 

    //- Yield of biomass (Y)
    yieldBiomass yB_;

    //- Yield of product (f)
    yieldProduct yP_;

    //- Carbon content
    carbonContent CC_;	

    //- Inhibition factor for h2 and nh3
    inhibitionParaI KI_;

    //- Inhibition factor
    inhibitionParaS KS_;

	//- Acid-base kinematic para [M/day]
    dimensionedScalar kAB_;

    //- Acid-base equilibrium para [M]
    admAcidBaseEq Ka_;

    //- Henry's Law coeffecients
    gasHenryCoe KH_;

    // Gas transfer rate coefficients [d^-1]
    dimensionedScalar kLa_;  

    // pH Limits
    limPH pHL_;

    //- Nitrogen content
    nitrogenContent NC_;

    //- Initial values
    List<scalar> Gini_;   // gas
    List<scalar> Mini_;   // medians
    scalar Pini_;         // initial for protons
    List<scalar> Eini_;   // electrolytes
    List<scalar> Iini_;   // ions

    //- DEBUG benchmark inflow
    // List<scalar> INFLOW_;
    List<dimensionedScalar> INFLOW_;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

    int printErrMessage(word runMode);

    dimensionedScalar defineTop(word runMode);

    decayRate defineRC(word runMode);

    yieldBiomass defineYields(word runMode);

    inhibitionParaI defineKI(word runMode);

    inhibitionParaS defineKS(word runMode);

    void defineSTOI();

    void defineInitialState(word runMode);

    // DEBUG
    void defineINFLOW(word runMode);

public:

    //- Biochemical rate coefficients (nu_ij)
    List<List<dimensionedScalar>> STOI;

    //- Some acied base stuff
    List<List<dimensionedScalar>> abDAE;

    //- Constructor

    admPara() : 
        errMessage
        (
            printErrMessage("Meso")
        ),
        ds_
        (
            dimMass/dimVolume
        ),
        Tbase_
        (
            298.15
        ),
        Top_
        (
            defineTop("Meso")
        ),
        kDec_
        (
            defineRC("Meso")
        ),
        yB_
        (
            defineYields("Meso")
        ),
        yP_
        (   
            0.10, // si_xc
            0.20, // 0.25, // xi_xc  <<< Rosen et al.
            0.20, // ch_xc
            0.20, // pr_xc
            0.30, // 0.25, // li_xc  <<< Rosen et al.
            0.95, // fa_li
            0.19, // h2_su
            0.13, // bu_su
            0.27, // pro_su
            0.41, // ac_su
            0.06, // h2_aa
            0.23, // va_aa
            0.26, // bu_aa
            0.05, // pro_aa
            0.40  // ac_aa
        ),
        CC_
        (
            0.2786, // xc
            0.03,   // si
            0.0313, // ch
            0.03,   // pr
            0.022,  // li
            0.03,   // xi
            0.0313, // su
            0.03,   // aa
            0.0217, // fa
            0.025,  // bu
            0.0268, // pro
            0.0313, // ac
            0.0313, // bac
            0.024,  // va
            0.0156  // ch4
        ),
        KI_
        (
            defineKI("Meso")
        ),
        KS_
        (
            defineKS("Meso")
        ),
        kAB_
        (
            dimMoles/dimTime, 1e8
        ),
        Ka_
        (
            1.380e-5,  // va
            1.514e-5,  // bu
            1.318e-5,  // pro
            1.738e-5,  // ac
            4.467e-7,  // co2
            5.623e-10, // IN
            1e-14      // W
        ),
        KH_
        (
            7.8e-04,    // h2
            0.0014,     // ch4
            0.035,      // co2
            0.0313      // h2o
        ),
        kLa_
        (
            dimless/dimTime, 200.0
        ),
        pHL_
        (
            5.5, 4.0, 7.0, 6.0, 6.0, 5.0
        ),
        NC_
        (
            0.0376/14.0, // xc
            0.06/14.0,   // I
            0.007,       // aa
            0.08/14.0    // bac
        )
    {
        defineInitialState("Meso");
        defineSTOI();
        defineINFLOW("Meso");
    };

    admPara
    (
        word runMode
    );

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void setParaDim(dimensionSet ds);

    inline scalar DTOS() const;

    inline const dimensionedScalar Tbase() const;

    inline const dimensionedScalar T() const;

    inline const decayRate kDec() const;

    inline const yieldBiomass yB() const;

    inline const yieldProduct yP() const;

    inline const carbonContent CC() const;

    inline const inhibitionParaI KI() const;

    inline const inhibitionParaS KS() const;

    inline const dimensionedScalar kAB() const;

    inline const admAcidBaseEq Ka() const;

    inline const gasHenryCoe KH() const;

    inline const dimensionedScalar kLa() const;

    inline const limPH pHL() const;

    inline const nitrogenContent NC() const;

    inline scalar Gini(label idx) const;

    inline scalar Mini(label idx) const;

    inline scalar Pini() const;

    inline scalar Eini(label idx) const;

    inline scalar Iini(label idx) const;

    // inline scalar INFLOW(label idx) const;

    inline dimensionedScalar INFLOW(label idx) const
    {
        return INFLOW_[idx];
    };
    
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#include "ADMno1ParameterI.H"
// #ifdef NoRepository
//     #include "ADMno1Templates.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //