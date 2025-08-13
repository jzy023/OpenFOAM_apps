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

//- Inhibition calculations

volScalarField::Internal Foam::ADMno1::calcInhibition
(
    volScalarField Y, 
    dimensionedScalar denum
)
{
    return 1.0/ (1.0 + (Y / denum));
}

volScalarField::Internal Foam::ADMno1::dCalcInhibition
(
    volScalarField Y, 
    dimensionedScalar denom
)
{
    return (-1 / (denom * (1 + (Y/denom)) * (1 + (Y/denom))));
}

volScalarField::Internal Foam::ADMno1::calcInhibitionHP
(
    volScalarField::Internal Shp,
    const dimensionedScalar UL,
    const dimensionedScalar LL,
    const dimensionedScalar n
)
{
    dimensionedScalar Kph = pow(10, -0.5 * (UL + LL));
    Kph.dimensions().reset(Shp.dimensions());

    return pow(Kph, n) / (pow(Shp / para_.MTOm(), n) + pow(Kph, n));
}


//- Kinetic rate calculations

volScalarField::Internal Foam::ADMno1::calcRho
(
    const dimensionedScalar k,
    volScalarField X
)
{
    return k * X.internalField();
}

volScalarField::Internal Foam::ADMno1::calcRho
(
    const dimensionedScalar k, 
    volScalarField S,
    const dimensionedScalar K,
    volScalarField X,
    volScalarField::Internal I
)
{
    return k * (S.internalField() / (K + S.internalField())) * X.internalField() * I;
}

volScalarField::Internal Foam::ADMno1::calcRho
(
    const dimensionedScalar k, 
    volScalarField S1,
    const dimensionedScalar K,
    volScalarField X,
    volScalarField S2,
    volScalarField::Internal I
)
{
    return k * (S1.internalField() / (K + S1.internalField())) * X.internalField() * 
           (1.0 / (1.0 + (S2.internalField() / S1.internalField()))) * I;
}


//- Components source term calculations

volScalarField::Internal Foam::ADMno1::concPerComponent
(
    label j,
    PtrList<volScalarField::Internal> KRPtrs
)
{
    volScalarField::Internal dY
    (
        IOobject
        (
            "dY",
            KRPtrs[0].mesh().time().timeName(),
            KRPtrs[0].mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        KRPtrs[0].mesh(),
        dimensionedScalar
        (
           "dY_Default", 
            KRPtrs[0].dimensions(),
            Zero
        )
    );

    // TODO: room for optimization (dont loop through all 19 elements since a lot of them are 0s)
    if (j == 9 || j == 10) // SIC or SIN
    {
        for (int i = 0; i < 12; i++) 
        {
            // dY += para_.DTOS() * KRPtrs[i] * para_.STOI[i][j];
            dY += para_.DTOS() * KRPtrs[i] * para_.STOI[i][j] * para_.MTOm(); 
        }
        for (int i = 12; i < 19; i++) 
        {
            // dY += para_.DTOS() * KRPtrs[i] * para_.STOI[12][j];
            dY += para_.DTOS() * KRPtrs[i] * para_.STOI[12][j] * para_.MTOm(); 
        }
    }
    else
    {
        for (int i = 0; i < 19; i++) 
        {
            dY += para_.DTOS() * KRPtrs[i] * para_.STOI[i][j]; 
        }
    }

    return dY;
}


//- Acid-base calculations

volScalarField::Internal Foam::ADMno1::fSion
(
    const dimensionedScalar Kax,
    const volScalarField::Internal& Sx,
    volScalarField::Internal Shp
)
{
    return Kax * Sx / (Kax + Shp);
}

volScalarField::Internal Foam::ADMno1::fSion
(
    const volScalarField::Internal& Kax,
    const volScalarField::Internal& Sx,
    volScalarField::Internal Shp
)
{
    return Kax * Sx / (Kax + Shp);
}

volScalarField::Internal Foam::ADMno1::dfSion
(
    const dimensionedScalar Kax,
    const volScalarField::Internal& Sx,
    volScalarField::Internal Shp
)
{
    return - Kax * Sx / ((Kax + Shp) * (Kax + Shp));
}

volScalarField::Internal Foam::ADMno1::dfSion
(
    const volScalarField::Internal& Kax,
    const volScalarField::Internal& Sx,
    volScalarField::Internal Shp
)
{
    return - Kax * Sx / ((Kax + Shp) * (Kax + Shp));
}


// ************************************************************************* //