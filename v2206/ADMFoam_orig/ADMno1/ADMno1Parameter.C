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
    C file for defining AMDno1 parameter class functions

\*---------------------------------------------------------------------------*/

#include "ADMno1Parameter.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

admPara::admPara
(
    word runMode
)
:
    errMessage
    (
        printErrMessage(runMode)
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
        defineTop(runMode)
    ),
    kDec_
    (
        defineRC(runMode)
    ),
    yB_
    (
        defineYields(runMode)
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
        // 0.2786, // xc
        0.02786, // xc
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
        defineKI(runMode)
    ),
    KS_
    (
        defineKS(runMode)
    ),
    kAB_
    (
        dimMoles/dimTime, 1e8
    ),
    Ka_
    (
        std::pow(10, -4.86),  // va
        std::pow(10, -4.82),  // bu
        std::pow(10, -4.88),  // pro
        std::pow(10, -4.76),  // ac
        std::pow(10, -6.35),  // co2
        std::pow(10, -9.25),  // IN
        std::pow(10, -14)     // W
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
    //  ULaa,LLaa,ULac,LLac,ULh2,LLh2
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
    // DEBUG
    defineINFLOW(runMode);
    Info<< "Ka_va:\t" << Ka_.va.value()
         << "\nKa_bu:\t" << Ka_.bu.value()
         << "\nKa_pro:\t" << Ka_.pro.value()
         << "\nKa_ac:\t" << Ka_.ac.value()
         << "\nKa_co2:\t" << Ka_.co2.value()
         << "\nKa_IN:\t" << Ka_.IN.value()
         << "\nKa_W:\t" << Ka_.W.value() << endl;

    defineInitialState(runMode);
    defineSTOI();
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

int admPara::printErrMessage
(
    word runMode
)
{
    if (!namesOpMode.found(runMode)) 
    {
        std::cerr << nl << "--> FOAM FATAL IO ERROR:" << nl
                  << "Unknown operation mode type: " << runMode
                  << "\n\nValid operation mode types:\n";
        std::cerr << namesOpMode.size() << "(\n";
        forAll(namesOpMode, i)
        {
            std::cerr << namesOpMode[i] << nl;
        }
        std::cerr << ")\n";
        std::exit(1);
    }

    return 0;
}

dimensionedScalar admPara::defineTop
(
    word runMode
)
{
    if(runMode == "Thermo")
    {
        return dimensionedScalar(dimTemperature, 273.15+55);
    }
    else
    {
        return dimensionedScalar(dimTemperature, 273.15+35);
    }
}

decayRate admPara::defineRC
(
    word runMode
)
{
    if(runMode == "Meso")
    {
        return decayRate
        (
            0.4,    // dis
            0.25,   // hyd_ch
            0.2,    // hyd_pr
            0.1,    // hyd_li
            30.0,   // m_su
            50.0,   // m_aa
            6.0,    // m_fa
            20.0,   // m_c4
            13.0,   // m_pro
            8.0,    // m_ac
            35.0,   // m_h2
            0.02,   // dec_xsu
            0.02,   // dec_xaa
            0.02,   // dec_xfa
            0.02,   // dec_xc4
            0.02,   // dec_xpro
            0.02,   // dec_xac
            0.02    // dec_xh2
        );
    }
    else if(runMode == "MesoSolid")
    {
        return decayRate
        (
            0.5,    // dis
            10.0,   // hyd_ch
            10.0,   // hyd_pr
            10.0,   // hyd_li
            30.0,   // m_su
            50.0,   // m_aa
            6.0,    // m_fa
            20.0,   // m_c4
            13.0,   // m_pro
            8.0,    // m_ac
            35.0,   // m_h2
            0.02,   // dec_xsu
            0.02,   // dec_xaa
            0.02,   // dec_xfa
            0.02,   // dec_xc4
            0.02,   // dec_xpro
            0.02,   // dec_xac
            0.02    // dec_xh2
        );
    }
    else
    {
        return decayRate
        (
            1.0,    // dis
            10.0,   // hyd_ch
            10.0,   // hyd_pr
            10.0,   // hyd_li
            70.0,   // m_su
            70.0,   // m_aa
            10.0,   // m_fa
            30.0,   // m_c4
            20.0,   // m_pro
            16.0,   // m_ac
            35.0,   // m_h2
            0.04,   // dec_xsu
            0.04,   // dec_xaa
            0.04,   // dec_xfa
            0.04,   // dec_xc4
            0.04,   // dec_xpro
            0.04,   // dec_xac
            0.04    // dec_xh2
        );
    }
}

yieldBiomass admPara::defineYields
(
    word runMode
)
{
    if(runMode != "Thermo")
    {
        return yieldBiomass
        (
            0.10, // su
            0.08, // aa
            0.06, // fa
            0.06, // c4
            0.04, // pro
            0.05, // ac
            0.06  // h2
        );
    }
    else
    {
        return yieldBiomass
        (
            0.10, // su
            0.08, // aa
            0.06, // fa
            0.06, // c4
            0.05, // pro
            0.05, // ac
            0.06  // h2
        );
    }
}

inhibitionParaI admPara::defineKI
(
    word runMode
)
{
    if(runMode == "Meso")
    {
        return inhibitionParaI
        (
            ds_,
            5.0e-6, // h2fa
	        1.0e-5, // h2c4
	        3.5e-6, // h2pro
	        0.0018  // nh3
        );
    }
    else if(runMode == "MesoSolid")
    {   // TODO: check value!
        return inhibitionParaI
        (
            ds_,
            5.0e-6, // h2fa
	        1.0e-5, // h2c4
	        3.5e-6, // h2pro
	        0.0018  // nh3
        );
    }
    else
    {
        return inhibitionParaI
        (   // from Rosen ADM-BSM DIGESTERPAR[55]
            ds_,
            5.0e-6, // h2fa
	        3.0e-5, // h2c4
	        1e-5,   // h2pro
	        0.011   // nh3
        );
    }
}

inhibitionParaS admPara::defineKS
(
    word runMode
)
{
    if(runMode == "Meso")
    {
        return inhibitionParaS
        (
            ds_,
            1e-4,   // IN
            1e-4,   // nh3
            0.5,    // su
            0.3,    // aa
            0.4,    // fa
            0.3,    // c4
            0.3,    // pro
            0.15,   // ac
            2.5e-5  // h2
        );
    }
    else if(runMode == "MesoSolid")
    {   // TODO: check value!
        return inhibitionParaS
        (
            ds_,
            1e-4,   // IN
            1e-4,   // nh3
            0.5,    // su
            0.3,    // aa
            0.4,    // fa
            0.2,    // c4
            0.1,    // pro
            0.15,   // ac
            7e-6    // h2
        );
    }
    else
    {
        return inhibitionParaS
        (
            ds_,
            1e-4,   // IN
            1e-4,   // nh3
            1.0,    // su
            0.3,    // aa
            0.4,    // fa
            0.4,    // c4
            0.3,    // pro
            0.3,    // ac
            5e-5    // h2
        );
    }
}

void admPara::defineSTOI()
{
    STOI.resize(19);
    for (int i = 0; i < 19; i++)
	{
	    STOI[i].resize(24);
	}
	// rows for processes and columns for substrates (visco[j,i] from the table)
    STOI[1][0] = 1.0;
    STOI[2][1] = 1.0;
    STOI[3][0] = 1.0- yP_.fa_li;
    STOI[3][2] = yP_.fa_li;

    STOI[0][11] = yP_.si_xc;
    STOI[0][12] = -1.0;
    STOI[0][13] = yP_.ch_xc;
    STOI[0][14] = yP_.pr_xc;
    STOI[0][15] = yP_.li_xc;
    STOI[0][23] = yP_.xi_xc;

    STOI[4][0] = -1.0;
    STOI[4][4] = (1.0- yB_.su) * yP_.bu_su;
    STOI[4][5] = (1.0- yB_.su) * yP_.pro_su;
    STOI[4][6] = (1.0- yB_.su) * yP_.ac_su;
    STOI[4][7] = (1.0- yB_.su) * yP_.h2_su;

    STOI[5][1] = -1.0;
    STOI[5][3] = (1.0- yB_.aa) * yP_.va_aa;
    STOI[5][4] = (1.0- yB_.aa) * yP_.bu_aa;
    STOI[5][5] = (1.0- yB_.aa) * yP_.pro_aa;
    STOI[5][6] = (1.0- yB_.aa) * yP_.ac_aa;
    STOI[5][7] = (1.0- yB_.aa) * yP_.h2_aa;

    STOI[6][2] = -1.0;
    STOI[6][6] = (1.0- yB_.fa) * 0.7;
    STOI[6][7] = (1.0- yB_.fa) * 0.3;

    STOI[7][3] = -1.0;
    STOI[7][5] = (1.0- yB_.c4) * 0.54;
    STOI[7][6] = (1.0- yB_.c4) * 0.31;
    STOI[7][7] = (1.0- yB_.c4) * 0.15;

    STOI[8][4] = -1.0;
    STOI[8][6] = (1.0- yB_.c4) * 0.8;
    STOI[8][7] = (1.0- yB_.c4) * 0.2;

    STOI[9][5] = -1.0;
    STOI[9][6] = (1.0- yB_.pro) * 0.57;
    STOI[9][7] = (1.0- yB_.pro) * 0.43;

    STOI[10][6] = -1.0;
    STOI[10][8] = (1.0- yB_.ac);

    STOI[11][7] = -1.0;
    STOI[11][8] = (1.0- yB_.h2);

    STOI[1][13] = -1.0;
    STOI[2][14] = -1.0;
    STOI[3][15] = -1.0;

    STOI[4][16] = yB_.su;
    STOI[5][17] = yB_.aa;
    STOI[6][18] = yB_.fa;
    STOI[7][19] = yB_.c4;

    STOI[8][19] = yB_.c4;
    STOI[9][20] = yB_.pro;
    STOI[10][21] = yB_.ac;
    STOI[11][22] = yB_.h2;

    for (int i = 12; i < 19; i++)
	{
	    STOI[i][12] = 1.0;
	    STOI[i][i + 4] = -1.0;
	}

	// Carbon content with Rosen paper
    STOI[0][9] = -(STOI[0][11] * CC_.si + STOI[0][12] * CC_.xc + STOI[0][13] * CC_.ch +
				   STOI[0][14] * CC_.pr + STOI[0][15] * CC_.li + STOI[0][23] * CC_.xi); // <<< added by Rosen et al.
    STOI[1][9] = -(STOI[1][0] * CC_.su +  STOI[1][13] * CC_.ch); // <<< added by Rosen et al.
    STOI[2][9] = -(STOI[2][1] * CC_.aa + STOI[2][14] * CC_.pr); // <<< added by Rosen et al.
    STOI[3][9] = -(STOI[3][0] * CC_.su + STOI[3][2] * CC_.fa + STOI[3][15] * CC_.li); // <<< added by Rosen et al.
    STOI[4][9] = -(STOI[4][0] * CC_.su + STOI[4][4] * CC_.bu + STOI[4][5] * CC_.pro +
	               STOI[4][6] * CC_.ac + STOI[4][16] * CC_.bac);
    STOI[5][9] = -(STOI[5][1] * CC_.aa + STOI[5][3] * CC_.va + STOI[5][4] * CC_.bu +
				   STOI[5][5] * CC_.pro + STOI[5][6] * CC_.ac + STOI[5][17] * CC_.bac);
    STOI[6][9] = -(STOI[6][2] * CC_.fa + STOI[6][6] * CC_.ac + STOI[6][18] * CC_.bac); // <<< added by Rosen et al.
    STOI[7][9] = -(STOI[7][3] * CC_.va + STOI[7][5] * CC_.pro + STOI[7][6] * CC_.ac + STOI[7][19] * CC_.bac); // <<< added by Rosen et al.
    STOI[8][9] = -(STOI[8][4] * CC_.bu + STOI[8][6] * CC_.ac + STOI[8][19] * CC_.bac); // <<< added by Rosen et al.
    STOI[9][9] = -(STOI[9][5] * CC_.pro + STOI[9][6] * CC_.ac + STOI[9][20] * CC_.bac);
    STOI[10][9] = -(STOI[10][6] * CC_.ac + STOI[10][8] * CC_.ch4 + STOI[10][21] * CC_.bac);
    STOI[11][9] = -(STOI[11][8] * CC_.ch4 + STOI[11][22] * CC_.bac);
    STOI[12][9] = -(STOI[12][12] * CC_.xc +  STOI[12][16] * CC_.bac); // <<< added by Rosen et al.

    // SIN calculation according to Rosen paper
    STOI[0][10] = NC_.xc - STOI[0][23] * NC_.I - STOI[0][11] * NC_.I - STOI[0][14] * NC_.aa; // <<< added by Rosen et al.
    STOI[4][10] = -STOI[4][16] * NC_.bac;
    STOI[5][10] = NC_.aa - STOI[5][17] * NC_.bac;
    STOI[6][10] = -STOI[6][18] * NC_.bac;
    STOI[7][10] = -STOI[7][19] * NC_.bac;
    STOI[8][10] = -STOI[8][19] * NC_.bac;
    STOI[9][10] = -STOI[9][20] * NC_.bac;
    STOI[10][10] = -STOI[10][21] * NC_.bac;
    STOI[11][10] = -STOI[11][22] * NC_.bac;
    STOI[12][10] = NC_.bac - NC_.xc; // <<< added by Rosen et al.

    for (label i = 0; i < 19; i++)
    {
        for (label j = 0; j < 24; j++)
        {
            STOI[i][j].dimensions().reset(dimless);
        }
    }
}


void admPara::defineInitialState(word runMode)
{
    Gini_.resize(3);
    Mini_.resize(3);
    Eini_.resize(5);
    Iini_.resize(2);

    Iini_[0] = 0.0;      // Scat
    Iini_[1] = 0.0052;   // San

    if(runMode == "Meso")
    {// TODO: update values
        Gini_[0] = 1.103241005083344e-05; // G_h2
        Gini_[1] = 1.653498470386505;     // G_ch4
        Gini_[2] = 0.013540127797408;     // G_co2
        Mini_[0] = 0.00942;               // S_co2
        Mini_[1] = 0.001884013268717;     // S_nh3
        Pini_ = 5.456176966644033e-08;    // S_hP

        Eini_[0] = 0.012283973156161;     // S_vaN
        Eini_[1] = 0.013952735474510;     // S_buN
        Eini_[2] = 0.017511435081451;     // S_proN
        Eini_[3] = 0.089035207194772;     // S_acN
        Eini_[4] = 0.085680011345346;     // S_hco3N
    }
    else if(runMode == "MesoSolid")
    {
        // Gini_[0] = 1.103241005083344e-05; // G_h2
        // Gini_[1] = 1.653498470386505;     // G_ch4
        // Gini_[2] = 0.013540127797408;     // G_co2
        // Mini_[0] = 0.00942;               // S_co2
        // Mini_[1] = 0.001884013268717;     // S_nh3
        // Pini_ = 5.456176966644033e-08;    // S_hP

        // Eini_[0] = 0.012283973156161;     // S_vaN
        // Eini_[1] = 0.013952735474510;     // S_buN
        // Eini_[2] = 0.017511435081451;     // S_proN
        // Eini_[3] = 0.089035207194772;     // S_acN
        // Eini_[4] = 0.085680011345346;     // S_hco3N

        Gini_[0] = 1.1032e-5; // G_h2
        Gini_[1] = 1.6535;    // G_ch4
        Gini_[2] = 0.0135;    // G_co2
        Mini_[0] = 0.00942;   // S_co2
        Mini_[1] = 0.001884;  // S_nh3
        Pini_ = 5.4562e-8;    // S_hP

        Eini_[0] = 0.012284;  // S_vaN
        Eini_[1] = 0.013953;  // S_buN
        Eini_[2] = 0.017511;  // S_proN
        Eini_[3] = 0.089035;  // S_acN
        Eini_[4] = 0.08568;   // S_hco3N
    }
    else
    {// TODO: update values
        Gini_[0] = 1.103241005083344e-05; // G_h2
        Gini_[1] = 1.653498470386505;     // G_ch4
        Gini_[2] = 0.013540127797408;     // G_co2
        Mini_[0] = 0.00942;               // S_co2
        Mini_[1] = 0.001884013268717;     // S_nh3
        Pini_ = 5.456176966644033e-08;    // S_hP

        Eini_[0] = 0.012283973156161;     // S_vaN
        Eini_[1] = 0.013952735474510;     // S_buN
        Eini_[2] = 0.017511435081451;     // S_proN
        Eini_[3] = 0.089035207194772;     // S_acN
        Eini_[4] = 0.085680011345346;     // S_hco3N
    }
}


void admPara::defineINFLOW
(
    word runMode
)
{
    INFLOW_.resize(26);

    if(runMode == "Meso")
    {   // TODO: update values
        INFLOW_[0] = .0;     // Ssu
        INFLOW_[1] = .043879921101364;  // Saa
        INFLOW_[2] = .0;     // Sfa
        INFLOW_[3] = .0;     // Sva
        INFLOW_[4] = .0;     // Sbu
        INFLOW_[5] = .0;     // Spro
        INFLOW_[6] = .0;     // Sac
        INFLOW_[7] = .0;     // Sh2
        INFLOW_[8] = .0;     // Sch4
        INFLOW_[9] = .007932590852686;    // SIC
        INFLOW_[10] = .001972071773301;   // SIN
        INFLOW_[11] = .028066505735355;   // SI
        INFLOW_[12] = .0;     // Xc
        INFLOW_[13] = 3.723594855179344;  // Xch
        INFLOW_[14] = 15.923520940540483; // Xpr
        INFLOW_[15] = 8.046980172357708;  // Xli
        INFLOW_[16] = .0;    // Xsu
        INFLOW_[17] = .0;    // Xaa
        INFLOW_[18] = .0;    // Xfa
        INFLOW_[19] = .0;    // Xc4
        INFLOW_[20] = .0;    // Xpro
        INFLOW_[21] = .0;    // Xac
        INFLOW_[22] = .0;    // Xh2
        INFLOW_[23] = 17.010642217805245; // XI
        INFLOW_[24] = .0;    // Scat+
        INFLOW_[25] = 0.005210099433331;  // San-
    }
    else if(runMode == "MesoSolid")
    {
        INFLOW_[0]  = 0.0;      // Ssu
        INFLOW_[1]  = 0.0439;   // Saa
        INFLOW_[2]  = 0.0;      // Sfa
        INFLOW_[3]  = 0.0;      // Sva		
        INFLOW_[4]  = 0.0;      // Sbu		
        INFLOW_[5]  = 0.0;      // Spro
        INFLOW_[6]  = 0.0;      // Sac		
        INFLOW_[7]  = 0.0;      // Sh2		
        INFLOW_[8]  = 0.0;      // Sch4
        INFLOW_[9]  = 0.0079;   // SIC
        INFLOW_[10] = 0.0020;   // SIN     //kmol/m3
        INFLOW_[11] = 0.0281;   // SI
        INFLOW_[12] = 0.0;      // Xc         //2.0
        INFLOW_[13] = 3.7236;   // Xch
        INFLOW_[14] = 15.9235;  // Xpr
        INFLOW_[15] = 8.0470;   // Xli
        INFLOW_[16] = 0.0;      // Xsu
        INFLOW_[17] = 0.0;      // Xaa         
        INFLOW_[18] = 0.0;      // Xfa			
        INFLOW_[19] = 0.0;      // Xc4
        INFLOW_[20] = 0.0;      // Xpro         
        INFLOW_[21] = 0.0;      // Xac			
        INFLOW_[22] = 0.0;      // Xh2
        INFLOW_[23] = 17.0106;  // XI
        INFLOW_[24] = 0.0;      // Scat+
        INFLOW_[25] = 0.00520;  // San-
        // INFLOW_[0] = .0;     // Ssu
        // INFLOW_[1] = .043879921101364;  // Saa
        // INFLOW_[2] = .0;     // Sfa
        // INFLOW_[3] = .0;     // Sva
        // INFLOW_[4] = .0;     // Sbu
        // INFLOW_[5] = .0;     // Spro
        // INFLOW_[6] = .0;     // Sac
        // INFLOW_[7] = .0;     // Sh2
        // INFLOW_[8] = .0;     // Sch4
        // INFLOW_[9] = .007932590852686;    // SIC
        // INFLOW_[10] = .001972071773301;   // SIN
        // INFLOW_[11] = .028066505735355;   // SI
        // INFLOW_[12] = .0;     // Xc
        // INFLOW_[13] = 3.723594855179344;  // Xch
        // INFLOW_[14] = 15.923520940540483; // Xpr
        // INFLOW_[15] = 8.046980172357708;  // Xli
        // INFLOW_[16] = .0;    // Xsu
        // INFLOW_[17] = .0;    // Xaa
        // INFLOW_[18] = .0;    // Xfa
        // INFLOW_[19] = .0;    // Xc4
        // INFLOW_[20] = .0;    // Xpro
        // INFLOW_[21] = .0;    // Xac
        // INFLOW_[22] = .0;    // Xh2
        // INFLOW_[23] = 17.010642217805245; // XI
        // INFLOW_[24] = .0;    // Scat+
        // INFLOW_[25] = 0.00520;  // San-
    }
    else
    {   // TODO: update values
        INFLOW_[0] = .0;     // Ssu
        INFLOW_[1] = .043879921101364;  // Saa
        INFLOW_[2] = .0;     // Sfa
        INFLOW_[3] = .0;     // Sva
        INFLOW_[4] = .0;     // Sbu
        INFLOW_[5] = .0;     // Spro
        INFLOW_[6] = .0;     // Sac
        INFLOW_[7] = .0;     // Sh2
        INFLOW_[8] = .0;     // Sch4
        INFLOW_[9] = .007932590852686;    // SIC
        INFLOW_[10] = .001972071773301;   // SIN
        INFLOW_[11] = .028066505735355;   // SI
        INFLOW_[12] = .0;     // Xc
        INFLOW_[13] = 3.723594855179344;  // Xch
        INFLOW_[14] = 15.923520940540483; // Xpr
        INFLOW_[15] = 8.046980172357708;  // Xli
        INFLOW_[16] = .0;    // Xsu
        INFLOW_[17] = .0;    // Xaa
        INFLOW_[18] = .0;    // Xfa
        INFLOW_[19] = .0;    // Xc4
        INFLOW_[20] = .0;    // Xpro
        INFLOW_[21] = .0;    // Xac
        INFLOW_[22] = .0;    // Xh2
        INFLOW_[23] = 17.010642217805245; // XI
        INFLOW_[24] = .0;    // Scat+
        INFLOW_[25] = 0.005210099433331;  // San-
    }
};

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void admPara::setParaDim
(
    dimensionSet ds
)
{
    ds_.reset(ds);

    KI_.setDimension(ds);

    KS_.setDimension(ds);

    Ka_.setDimension(ds);

    for(label i = 0; i < 26; i++)
    {
        INFLOW_[i].dimensions().reset(ds);
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //