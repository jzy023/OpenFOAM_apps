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

// #include "fvCFD.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Inhibition calculations

volScalarField::Internal calcInhibition
(
    volScalarField Y, 
    dimensionedScalar denom
);

volScalarField::Internal dCalcInhibition
(
    volScalarField Y, 
    dimensionedScalar denom
);

volScalarField::Internal calcInhibitionHP
(
    volScalarField::Internal Shp,
    const dimensionedScalar UL,
    const dimensionedScalar LL,
    const dimensionedScalar n
);

//- Kinetic rate calculations

volScalarField::Internal calcRho
(
    const dimensionedScalar k, 
    volScalarField X
);

volScalarField::Internal calcRho
(
    const dimensionedScalar k, 
    volScalarField S,
    const dimensionedScalar K,
    volScalarField X,
    volScalarField::Internal I
);

volScalarField::Internal calcRho
(
    const dimensionedScalar k, 
    volScalarField S1,
    const dimensionedScalar K,
    volScalarField X,
    volScalarField S2,
    volScalarField::Internal I
);

//- Components source term calculations

volScalarField::Internal concPerComponent
(
    label j,
    PtrList<volScalarField::Internal> KRPtrs
);


//- Acid-base calculations

volScalarField::Internal fSion
(
    const dimensionedScalar Kax,
    const volScalarField::Internal& Sx,
    volScalarField::Internal Shp
);

volScalarField::Internal fSion
(
    const volScalarField::Internal& Kax,
    const volScalarField::Internal& Sx,
    volScalarField::Internal Shp
);

volScalarField::Internal dfSion
(
    const dimensionedScalar Kax,
    const volScalarField::Internal& Sx,
    volScalarField::Internal Shp
);

volScalarField::Internal dfSion
(
    const volScalarField::Internal& Kax,
    const volScalarField::Internal& Sx,
    volScalarField::Internal Shp
);


// ************************************************************************* //