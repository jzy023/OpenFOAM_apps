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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Sh2 calculations

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
    volScalarField::Internal GRSh2Temp = para_.DTOS() * para_.kLa() 
                                       * (Sh2Temp.internalField() - R_ * TopDummy_.internalField() * GPtrs_[0].internalField() * KHh2_);

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

    dKRPtrs_temp[11] = para_.kDec().m_h2 * YPtrs_[22].internalField() * IPtrs_[2] * IPtrs_[3] * para_.KS().h2
                     / ((para_.KS().h2 + Sh2Temp.internalField()) * (para_.KS().h2 + Sh2Temp.internalField()));

    // volScalarField dConv(fvc::div(flux));
    dimensionedScalar dConv = - para_.DTOS() * (Qin_/Vliq_);
    dimensionedScalar dGRSh2Temp = para_.DTOS() * para_.kLa();

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

//- Other calculations

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
    for(label j = 0; j < 7; j++)
    {
        volScalarField::Internal inOutFlow(dYPtrs_[0]);
        inOutFlow = para_.DTOS() * (Qin_/Vliq_) * (para_.INFLOW(j) - YPtrs_[j]);
        dYPtrs_[j] = concPerComponent(j, KRPtrs_) + inOutFlow;
        // dYPtrs_[j] = concPerComponent(j, KRPtrs_);
    }

    for(label j = 8; j < YPtrs_.size(); j++)
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

// ************************************************************************* //