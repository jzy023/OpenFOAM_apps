/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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
    Templated abstract base class for volumen occupancy incompressible
    turbulence models.

SourceFiles
    DPMIncompressibleTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef typeDefs_H
#define typeDefs_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// namespace ADM
// {

/** =======================================================================
 * STOI Components
 * @brief	ADM-1 Table 6.2, raction first order prameters [d^-1]
 */
struct reactionRateConstant
{
	scalar dis{};	//[1/d]
	scalar hyd_ch{}; //[1/d]
	scalar hyd_pr{}; //[1/d]
	scalar hyd_li{}; //[1/d]

	scalar m_su{};   //[1/d]
	scalar m_aa{};   //[1/d]
	scalar m_fa{};   //[1/d]
	scalar m_c4{};   //[1/d]
	scalar m_pro{};  //[1/d]
	scalar m_ac{};   //[1/d]
	scalar m_h2{};   //[1/d]

	scalar dec_xsu{};    //[1/d]
	scalar dec_xaa{};    //[1/d]
	scalar dec_xfa{};    //[1/d]
	scalar dec_xc4{};    //[1/d]
	scalar dec_xpro{};   //[1/d]
	scalar dec_xac{};    //[1/d]
	scalar dec_xh2{};    //[1/d]
};

/** =======================================================================
 * @brief coefficients for gas mass transfer [d^-1]
 *
 */
struct gasTransferRate
{
	scalar h2{};
	scalar ch4{};
	scalar co2{};
};

/** =======================================================================
 * @brief Henry's Law coefficients [M bar^-1]
 *
 */
struct gasHenryCoe
{
	scalar h2{};
	scalar ch4{};
	scalar co2{};
    scalar h2o{};
};

/** =======================================================================
 * @brief	ADM-1 Table 6.2, default Yields of biomass (Y) [COD/COD]
 */
struct yieldBiomass
{
	scalar su{};  //[COD/COD]
	scalar aa{};  //[COD/COD]
	scalar fa{};  //[COD/COD]
	scalar c4{};  //[COD/COD]
	scalar pro{}; //[COD/COD]
	scalar ac{};  //[COD/COD]
	scalar h2{};  //[COD/COD]
};

/** =======================================================================
 * @brief	ADM-1 Table 6.1, default Yields of product (f) [kgCOD/kgCOD]
 */
struct yieldProduct
{
	scalar si_xc{};  //[kgCOD/kgCOD]
	scalar xi_xc{};  //[kgCOD/kgCOD]
	scalar ch_xc{};  //[kgCOD/kgCOD]
	scalar pr_xc{};  //[kgCOD/kgCOD]
	scalar li_xc{};  //[kgCOD/kgCOD]
	scalar fa_li{};  //[kgCOD/kgCOD]
	scalar h2_su{};  //[kgCOD/kgCOD]
	scalar bu_su{};  //[kgCOD/kgCOD]
	scalar pro_su{}; //[kgCOD/kgCOD]
	scalar ac_su{};  //[kgCOD/kgCOD]
	scalar h2_aa{};  //[kgCOD/kgCOD]
	scalar va_aa{};  //[kgCOD/kgCOD]
	scalar bu_aa{};  //[kgCOD/kgCOD]
	scalar pro_aa{}; //[kgCOD/kgCOD]
	scalar ac_aa{};  //[kgCOD/kgCOD]
};

/** =======================================================================
 * @brief	default Carbon contents of substrates
 * 			_bac and _biom are interchangeable
 * 			Kmole/KgCOD
 */
struct carbonContent
{
	scalar xc{};  //[kmol/kgCOD]
	scalar si{};  //[kmol/kgCOD]
	scalar ch{};  //[kmol/kgCOD]
	scalar pr{};  //[kmol/kgCOD]
	scalar li{};  //[kmol/kgCOD]
	scalar xi{};  //[kmol/kgCOD]
	scalar su{};  //[kmol/kgCOD]
	scalar aa{};  //[kmol/kgCOD]
	scalar fa{};  //[kmol/kgCOD]
	scalar bu{};  //[kmol/kgCOD]
	scalar pro{}; //[kmol/kgCOD]
	scalar ac{};  //[kmol/kgCOD]
	scalar bac{}; //[kmol/kgCOD]
	scalar va{};  //[kmol/kgCOD]
	scalar ch4{}; //[kmol/kgCOD]
};

/** =======================================================================
 * Inhibition calcilations
 * @brief	default Inhibition parameters for (KI)
 */
struct inhibitionParaI
{
	scalar h2fa{};   // [KgCOD/m3] Fatty Acids
	scalar h2c4{};   // [KgCOD/m3] Valerate and Butyrate
	scalar h2pro{};  // [KgCOD/m3] Propionate
	scalar nh3{};    // [M=Kmole/m3] Ammonium
};

/** =======================================================================
 * @brief	default Inhibition parameters for [Ks]
 */
struct inhibitionParaS
{
	scalar IN{};  //[kmol/m3]
	scalar nh3{}; //[kmol/m3]
	scalar su{};  //[kgCOD/m3]
	scalar aa{};  //[kgCOD/m3]
	scalar fa{};  //[kgCOD/m3]
	scalar c4{};  //[kgCOD/m3]
	scalar pro{}; //[kgCOD/m3]
	scalar ac{};  //[kgCOD/m3]
	scalar h2{};  //[kgCOD/m3]
};

/** =======================================================================
 * @brief
 * Inhibition parameters, I
 */
struct inhibitionCoefI
{
	scalar phaa{};
	scalar phac{};
	scalar phh2{};
	scalar h2fa{};
	scalar h2c4{};
	scalar h2pro{};
	scalar nh3{};
	scalar IN{};
};

/** =======================================================================
 * @brief structure to keep the concetrations of Dynamic state variables
 */
struct admTankConc
{
	scalar S_su{};	 //[kg/m3]
	scalar S_aa{};	 //[kg/m3]
	scalar S_fa{};	 //[kg/m3]
	scalar S_va{};	 //[kg/m3]
	scalar S_bu{};	 //[kg/m3]
	scalar S_pro{};   //[kg/m3]
	scalar S_ac{};	 //[kg/m3]
	scalar S_h2{};	 //[kg/m3]
	scalar S_ch4{};   //[kg/m3]
	scalar S_IC{};	 //[kmol/m3]
	scalar S_IN{};	 //[Kmol/m3]
	scalar S_I{};	 //[kg/m3]

	scalar X_c{};	 //[kg/m3]
	scalar X_ch{};	 //[kg/m3]
	scalar X_pr{};	 //[kg/m3]
	scalar X_li{};	 //[kg/m3]
	scalar X_su{};	 //[kg/m3]
	scalar X_aa{};	 //[kg/m3]
	scalar X_fa{};	 //[kg/m3]
	scalar X_c4{};	 //[kg/m3]
	scalar X_pro{};   //[kg/m3]
	scalar X_ac{};	 //[kg/m3]
	scalar X_h2{};	 //[kg/m3]
	scalar X_I{};	 //[kg/m3]

	//==============================================================
	scalar S_cat{};   //[kmol/m3]
	scalar S_an{};	 //[kmol/m3]
	scalar S_hP{};	 //[kmol/m3]
	scalar S_vaN{};   //[kg/m3]
	scalar S_buN{};   //[kg/m3]
	scalar S_proN{};  //[kg/m3]
	scalar S_acN{};   //[kg/m3]
	scalar S_hco3N{}; //[kmol/m3]
	scalar S_nh4{};   //[kmol/m3]
	scalar S_nh3{};   //[kmol/m3]
	scalar S_ohN{};   //[kg/m3]
	// There is nh4+ as well, dont forget that during calculations

	scalar S_co2{};

	scalar G_h2{};	 // gas phase component [kg/m3]
	scalar G_ch4{};   // gas phase component [kg/m3]
	scalar G_co2{};   // gas phase component [kg/m3]
};

/** =======================================================================
 * @brief structure to keep the Diffusion concetrations of Dynamic state variables
 * should replace inflow
 */
struct admConcDiff
{
	scalar S_su{};	 //[kg/m3]
	scalar S_aa{};	 //[kg/m3]
	scalar S_fa{};	 //[kg/m3]
	scalar S_va{};	 //[kg/m3]
	scalar S_bu{};	 //[kg/m3]
	scalar S_pro{};   //[kg/m3]
	scalar S_ac{};	 //[kg/m3]
	scalar S_h2{};	 //[kg/m3]
	scalar S_ch4{};   //[kg/m3]
	scalar S_IC{};	 //[kmol/m3]
	scalar S_IN{};	 //[Kmol/m3]
	scalar S_I{};	 //[kg/m3]

    scalar G_h2{};    // gas phase component [kg/m3]
    scalar G_ch4{};   // gas phase component [kg/m3]
    scalar G_co2{};   // gas phase component [kg/m3]
};

/** =======================================================================
 * @brief structure to keep the inflow concetrations of Dynamic state variables
 */
struct admInflowConc
{
	scalar S_su{};	 //[kg/m3]
	scalar S_aa{};	 //[kg/m3]
	scalar S_fa{};	 //[kg/m3]
	scalar S_va{};	 //[kg/m3]
	scalar S_bu{};	 //[kg/m3]
	scalar S_pro{};   //[kg/m3]
	scalar S_ac{};	 //[kg/m3]
	scalar S_h2{};	 //[kg/m3]
	scalar S_ch4{};   //[kg/m3]
	scalar S_IC{};	 //[kmol/m3]
	scalar S_IN{};	 //[Kmol/m3]
	scalar S_I{};	 //[kg/m3]

	scalar X_c{};	 //[kg/m3]
	scalar X_ch{};	 //[kg/m3]
	scalar X_pr{};	 //[kg/m3]
	scalar X_li{};	 //[kg/m3]
	scalar X_su{};	 //[kg/m3]
	scalar X_aa{};	 //[kg/m3]
	scalar X_fa{};	 //[kg/m3]
	scalar X_c4{};	 //[kg/m3]
	scalar X_pro{};   //[kg/m3]
	scalar X_ac{};	 //[kg/m3]
	scalar X_h2{};	 //[kg/m3]
	scalar X_I{};	 //[kg/m3]

	//==============================================================
	scalar S_cat{};   //[kmol/m3]
	scalar S_an{};	 //[kmol/m3]

	// scalar S_hP;		//[kg/m3]
	// scalar S_vaN;		//[kg/m3]
	// scalar S_buN;		//[kg/m3]
	// scalar S_proN;		//[kg/m3]
	// scalar S_acN;		//[kg/m3]
	// scalar S_hco3N;	//[kmol/m3]
	// scalar S_nh4;		//[kmol/m3]
	// scalar S_nh3;		//[kmol/m3]
};

/** =======================================================================
 * @brief transfer rate of gas phase [ d^-1]
 *
 */
struct rhoGas
{
	scalar h2{};     // [kg/d]
	scalar ch4{};    // [kg/d]
	scalar co2{};    // [kg/d]
};

/** =======================================================================
 * @brief amount of gas exiting the cstr/particle
 *
 */
struct admGasExit
{
    scalar h2{};     // [kg/m3]
    scalar ch4{};    // [kg/m3]
    scalar co2{};    // [kg/m3]
};

/** =======================================================================
 * @brief structure to keep the Acide base Kinematic Para [kAB][/M/day]
 *
 */
struct admAcidBaseKinetic
{
	scalar va{};  //[m3/kmol/day]
	scalar bu{};  //[m3/kmol/day]
	scalar pro{}; //[m3/kmol/day]
	scalar ac{};  //[m3/kmol/day]
	scalar co2{}; //[m3/kmol/day]
	scalar IN{};  //[m3/kmol/day]
};

/** =======================================================================
 * @brief structure to keep the Acide base Equilibrium Para [Ka][M]
 */
struct admAcidBaseEq
{
	scalar va{};  //[M]
	scalar bu{};  //[M]
	scalar pro{}; //[M]
	scalar ac{};  //[M]
	scalar co2{}; //[M]
	scalar IN{};  //[M]
	scalar W{};   //[M]
};

// } // End namespace ADM

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //