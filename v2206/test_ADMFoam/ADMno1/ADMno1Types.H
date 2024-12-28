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
    ADM::typeDefs

Description
    Data structure definition for ADMno1 parameters .


\*---------------------------------------------------------------------------*/

#ifndef ADMno1Types_H
#define ADMno1Types_H

#include <iomanip>
#include "IOdictionary.H"
#include "typeInfo.H"
#include "fvCFD.H"
#include "autoPtr.H"
#include "ODESystem.H"
#include "volFields.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TODO: try making things std::unordered_map to clean up using for loops



//- STOI Components
// ADM-1 Table 6.2, raction first order prameters [d^-1]
struct decayRate
{
    dimensionedScalar dis;	   //[1/d]
    dimensionedScalar hyd_ch; //[1/d]
    dimensionedScalar hyd_pr; //[1/d]
    dimensionedScalar hyd_li; //[1/d]
    dimensionedScalar m_su;   //[1/d]
    dimensionedScalar m_aa;   //[1/d]
    dimensionedScalar m_fa;   //[1/d]
    dimensionedScalar m_c4;   //[1/d]
    dimensionedScalar m_pro;  //[1/d]
    dimensionedScalar m_ac;   //[1/d]
    dimensionedScalar m_h2;   //[1/d]
    dimensionedScalar dec_xsu;    //[1/d]
    dimensionedScalar dec_xaa;    //[1/d]
    dimensionedScalar dec_xfa;    //[1/d]
    dimensionedScalar dec_xc4;    //[1/d]
    dimensionedScalar dec_xpro;   //[1/d]
    dimensionedScalar dec_xac;    //[1/d]
    dimensionedScalar dec_xh2;    //[1/d]

    decayRate
    (
        scalar dis_,
        scalar hyd_ch_,
        scalar hyd_pr_,
        scalar hyd_li_,
        scalar m_su_,
        scalar m_aa_,
        scalar m_fa_,
        scalar m_c4_,
        scalar m_pro_,
        scalar m_ac_,
        scalar m_h2_,
        scalar dec_xsu_,
        scalar dec_xaa_,
        scalar dec_xfa_,
        scalar dec_xc4_,
        scalar dec_xpro_,
        scalar dec_xac_,
        scalar dec_xh2_
    )
    :
        dis( dimless/dimTime, dis_),
        hyd_ch( dimless/dimTime, hyd_ch_),
        hyd_pr( dimless/dimTime, hyd_pr_),
        hyd_li( dimless/dimTime, hyd_li_),
        m_su( dimless/dimTime, m_su_),
        m_aa( dimless/dimTime, m_aa_),
        m_fa( dimless/dimTime, m_fa_),
        m_c4( dimless/dimTime, m_c4_),
        m_pro( dimless/dimTime, m_pro_),
        m_ac( dimless/dimTime, m_ac_),
        m_h2( dimless/dimTime, m_h2_),
        dec_xsu( dimless/dimTime, dec_xsu_),
        dec_xaa( dimless/dimTime, dec_xaa_),
        dec_xfa( dimless/dimTime, dec_xfa_),
        dec_xc4( dimless/dimTime, dec_xc4_),
        dec_xpro( dimless/dimTime, dec_xpro_),
        dec_xac( dimless/dimTime, dec_xac_),
        dec_xh2( dimless/dimTime, dec_xh2_)
    {};
};

//- Henry's Law coefficients [M bar^-1]
// TODO: check dimension
struct gasHenryCoe
{
    dimensionedScalar h2;
    dimensionedScalar ch4;
    dimensionedScalar co2;
    dimensionedScalar h2o; // TODO: not in use

    gasHenryCoe
    (
        scalar h2_,
        scalar ch4_,
        scalar co2_,
        scalar h2o_
    )
    :
        h2( dimless, h2_ ),
        ch4( dimless, ch4_ ),
        co2( dimless, co2_ ),
        h2o( dimless, h2o_ )
    {};

};

//- ADM-1 Table 6.2, default Yields of biomass (Y) [COD/COD]
struct yieldBiomass
{
    dimensionedScalar su;  //[COD/COD]
    dimensionedScalar aa;  //[COD/COD]
    dimensionedScalar fa;  //[COD/COD]
    dimensionedScalar c4;  //[COD/COD]
    dimensionedScalar pro; //[COD/COD]
    dimensionedScalar ac;  //[COD/COD]
    dimensionedScalar h2;  //[COD/COD]

    yieldBiomass
    (
        scalar su_,
        scalar aa_,
        scalar fa_,
        scalar c4_,
        scalar pro_,
        scalar ac_,
        scalar h2_
    )
    :
        su( dimless, su_ ),
        aa( dimless, aa_ ),
        fa( dimless, fa_ ),
        c4( dimless, c4_ ),
        pro( dimless, pro_ ),
        ac( dimless, ac_ ),
        h2( dimless, h2_ )
    {};   
};

//- ADM-1 Table 6.1, default Yields of product (f) [kgCOD/kgCOD]
struct yieldProduct
{
    dimensionedScalar si_xc;  //[kgCOD/kgCOD]
    dimensionedScalar xi_xc;  //[kgCOD/kgCOD]
    dimensionedScalar ch_xc;  //[kgCOD/kgCOD]
    dimensionedScalar pr_xc;  //[kgCOD/kgCOD]
    dimensionedScalar li_xc;  //[kgCOD/kgCOD]
    dimensionedScalar fa_li;  //[kgCOD/kgCOD]
    dimensionedScalar h2_su;  //[kgCOD/kgCOD]
    dimensionedScalar bu_su;  //[kgCOD/kgCOD]
    dimensionedScalar pro_su; //[kgCOD/kgCOD]
    dimensionedScalar ac_su;  //[kgCOD/kgCOD]
    dimensionedScalar h2_aa;  //[kgCOD/kgCOD]
    dimensionedScalar va_aa;  //[kgCOD/kgCOD]
    dimensionedScalar bu_aa;  //[kgCOD/kgCOD]
    dimensionedScalar pro_aa; //[kgCOD/kgCOD]
    dimensionedScalar ac_aa;  //[kgCOD/kgCOD]

    yieldProduct
    (
        scalar si_xc_,
        scalar xi_xc_,
        scalar ch_xc_,
        scalar pr_xc_,
        scalar li_xc_,
        scalar fa_li_,
        scalar h2_su_,
        scalar bu_su_,
        scalar pro_su_,
        scalar ac_su_,
        scalar h2_aa_,
        scalar va_aa_,
        scalar bu_aa_,
        scalar pro_aa_,
        scalar ac_aa_
    )
    :
        si_xc( dimless, si_xc_ ),
        xi_xc( dimless, xi_xc_ ),
        ch_xc( dimless, ch_xc_ ),
        pr_xc( dimless, pr_xc_ ),
        li_xc( dimless, li_xc_ ),
        fa_li( dimless, fa_li_ ),
        h2_su( dimless, h2_su_ ),
        bu_su( dimless, bu_su_ ),
        pro_su( dimless, pro_su_ ),
        ac_su( dimless, ac_su_ ),
        h2_aa( dimless, h2_aa_ ),
        va_aa( dimless, va_aa_ ),
        bu_aa( dimless, bu_aa_ ),
        pro_aa( dimless, pro_aa_ ),
        ac_aa( dimless, ac_aa_ )
    {};
};

//- Default Carbon contents of substrates
//  Kmole/KgCOD
//  (_bac and _biom in literature are interchangeable )
struct carbonContent
{
    dimensionedScalar xc;  //[kmol/kgCOD]
    dimensionedScalar si;  //[kmol/kgCOD]
    dimensionedScalar ch;  //[kmol/kgCOD]
    dimensionedScalar pr;  //[kmol/kgCOD]
    dimensionedScalar li;  //[kmol/kgCOD]
    dimensionedScalar xi;  //[kmol/kgCOD]
    dimensionedScalar su;  //[kmol/kgCOD]
    dimensionedScalar aa;  //[kmol/kgCOD]
    dimensionedScalar fa;  //[kmol/kgCOD]
    dimensionedScalar bu;  //[kmol/kgCOD]
    dimensionedScalar pro; //[kmol/kgCOD]
    dimensionedScalar ac;  //[kmol/kgCOD]
    dimensionedScalar bac; //[kmol/kgCOD]
    dimensionedScalar va;  //[kmol/kgCOD]
    dimensionedScalar ch4; //[kmol/kgCOD]

    carbonContent
    (
        scalar xc_,
        scalar si_,
        scalar ch_,
        scalar pr_,
        scalar li_,
        scalar xi_,
        scalar su_,
        scalar aa_,
        scalar fa_,
        scalar bu_,
        scalar pro_,
        scalar ac_,
        scalar bac_,
        scalar va_,
        scalar ch4_
    )
    :
        xc( dimless, xc_ ),
        si( dimless, si_ ),
        ch( dimless, ch_ ),
        pr( dimless, pr_ ),
        li( dimless, li_ ),
        xi( dimless, xi_ ),
        su( dimless, su_ ),
        aa( dimless, aa_ ),
        fa( dimless, fa_ ),
        bu( dimless, bu_ ),
        pro( dimless, pro_ ),
        ac( dimless, ac_ ),
        bac( dimless, bac_ ),
        va( dimless, va_ ),
        ch4( dimless, ch4_ )
    {};
};

//- Default Inhibition parameters for (KI)
struct inhibitionParaI
{
    dimensionSet ds;

    dimensionedScalar h2fa;     // [KgCOD/m3] 
    dimensionedScalar h2c4;     // [KgCOD/m3] 
    dimensionedScalar h2pro;    // [KgCOD/m3] 
    dimensionedScalar nh3;      // [M=Kmole/m3]

    inhibitionParaI
    (
        dimensionSet ds_,
        scalar h2fa_,
        scalar h2c4_,
        scalar h2pro_,
        scalar nh3_
    )
    :
        h2fa( ds_, h2fa_ ),
        h2c4( ds_, h2c4_ ),
        h2pro( ds_, h2pro_ ),
        nh3( ds_, nh3_ )
    {};

    void setDimension
    (
        dimensionSet ds
    )
    {
        h2fa.dimensions().reset(ds);
        h2c4.dimensions().reset(ds);
        h2pro.dimensions().reset(ds);
        nh3.dimensions().reset(ds);
    }
};

//- Default Inhibition parameters for [Ks]
struct inhibitionParaS
{
    dimensionSet ds;

    dimensionedScalar IN;  //[kmol/m3]
    dimensionedScalar su;  //[kgCOD/m3]
    dimensionedScalar aa;  //[kgCOD/m3]
    dimensionedScalar fa;  //[kgCOD/m3]
    dimensionedScalar c4;  //[kgCOD/m3]
    dimensionedScalar pro; //[kgCOD/m3]
    dimensionedScalar ac;  //[kgCOD/m3]
    dimensionedScalar h2;  //[kgCOD/m3]

    inhibitionParaS
    (
        dimensionSet ds_,
        scalar IN_,
        scalar su_,
        scalar aa_,
        scalar fa_,
        scalar c4_,
        scalar pro_,
        scalar ac_,
        scalar h2_
    )
    :
        IN( ds_, IN_ ),
        su( ds_, su_ ),
        aa( ds_, aa_ ),
        fa( ds_, fa_ ),
        c4( ds_, c4_ ),
        pro( ds_, pro_ ),
        ac( ds_, ac_ ),
        h2( ds_, h2_ )
    {};

    void setDimension
    (
        dimensionSet ds
    )
    {
        IN.dimensions().reset(ds);
        su.dimensions().reset(ds);
        aa.dimensions().reset(ds);
        fa.dimensions().reset(ds);
        c4.dimensions().reset(ds);
        pro.dimensions().reset(ds);
        ac.dimensions().reset(ds);
        h2.dimensions().reset(ds);
    }
};

//- Structure to keep the Acid base Kinematic Para [kAB_][/M/day]
struct admAcidBaseKinetic
{
    dimensionedScalar va;  //[m3/kmol/day]
    dimensionedScalar bu;  //[m3/kmol/day]
    dimensionedScalar pro; //[m3/kmol/day]
    dimensionedScalar ac;  //[m3/kmol/day]
    dimensionedScalar co2; //[m3/kmol/day]
    dimensionedScalar IN;  //[m3/kmol/day]

    admAcidBaseKinetic
    (
        scalar va_,
        scalar bu_,
        scalar pro_,
        scalar ac_,
        scalar co2_,
        scalar IN_
    )
    :
        va( dimVolume/dimMass/dimTime, va_ ),
        bu( dimVolume/dimMass/dimTime, bu_ ),
        pro( dimVolume/dimMass/dimTime, pro_ ),
        ac( dimVolume/dimMass/dimTime, ac_ ),
        co2( dimVolume/dimMass/dimTime, co2_ ),
        IN( dimVolume/dimMass/dimTime, IN_ )
    {};
};

//- Structure to keep the Acid base Equilibrium Para [M]
struct admAcidBaseEq
{
    dimensionSet ds;

    dimensionedScalar va;  //[M]
    dimensionedScalar bu;  //[M]
    dimensionedScalar pro; //[M]
    dimensionedScalar ac;  //[M]
    dimensionedScalar co2; //[M]
    dimensionedScalar IN;  //[M]
    dimensionedScalar W;   //[M]

    admAcidBaseEq
    (
        scalar va_,
        scalar bu_,
        scalar pro_,
        scalar ac_,
        scalar co2_,
        scalar IN_,
        scalar W_
    )
    :
        va( dimMass, va_),
        bu( dimMass, bu_),
        pro( dimMass, pro_),
        ac( dimMass, ac_),
        co2( dimMass, co2_),
        IN( dimMass, IN_),
        W( dimMass, W_)
    {};

    void setDimension
    (
        dimensionSet ds
    )
    {
        va.dimensions().reset(ds);
        bu.dimensions().reset(ds);
        pro.dimensions().reset(ds);
        ac.dimensions().reset(ds);
        co2.dimensions().reset(ds);
        IN.dimensions().reset(ds);
        W.dimensions().reset(ds);
    }
};

//- pH Limits 
struct limPH
{
    dimensionedScalar ULaa;
    dimensionedScalar LLaa;
    dimensionedScalar ULac;
    dimensionedScalar LLac;
    dimensionedScalar ULh2;
    dimensionedScalar LLh2;

    limPH
    (
        scalar ULaa_,
        scalar LLaa_,
        scalar ULac_,
        scalar LLac_,
        scalar ULh2_,
        scalar LLh2_
    )
    :
        ULaa( dimless, ULaa_),
        LLaa( dimless, LLaa_),
        ULac( dimless, ULac_),
        LLac( dimless, LLac_),
        ULh2( dimless, ULh2_),
        LLh2( dimless, LLh2_)
    {};
};

//- Nitrogen content
// TODO: check dimensions
struct nitrogenContent
{
    dimensionedScalar xc;
    dimensionedScalar I;
    dimensionedScalar aa;
    dimensionedScalar bac; // biomass

    nitrogenContent
    (
        scalar xc_,
        scalar I_,
        scalar aa_,
        scalar bac_
    )
    :
        xc( dimless, xc_ ),
        I( dimless, I_ ),
        aa( dimless, aa_ ),
        bac( dimless, bac_)
    {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //