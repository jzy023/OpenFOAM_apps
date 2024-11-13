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
 
#include "ADMno1.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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

    volScalarField::Internal E = IOPtrs_[0].internalField() - IOPtrs_[1].internalField() + ShpTemp - SohN + (YPtrs_[10].internalField() - MPtrs_[1].internalField()) - 
                                 EPtrs_[4] - EPtrs_[3]/64.0 - EPtrs_[2]/112.0 - EPtrs_[1]/160.0 - EPtrs_[0]/208.0;

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

    return uniField - dSnh3 - dShco3N - dSacN/64.0 - 
           dSproN/112.0 - dSbuN/160.0 - dSvaN/208.0 - dSohN;
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
        max(mag(E.field())) > tol &&
        i < nIter
    );

    if( min(x.field()) < 0 )
    {
        x.field() = 0.0*x.field() + 1e-16;
    }

    Info<< "Newton-Raphson:\tSolving for Sh+" 
        << ", min Shp: " << min(x.field()) 
        << ", max Shp: " << max(x.field()) 
        << ", No Interations " << i << endl;

    // ShP
    ShP_ = x;
    pH_.field() = -log10(ShP_.field());

    // update Sco2: Sco2 = SIC - Shco3N
    MPtrs_[0].ref() = YPtrs_[9] - EPtrs_[4];

}


// ************************************************************************* //