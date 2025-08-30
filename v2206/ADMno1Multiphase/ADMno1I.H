/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |  
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

inline scalar 
Foam::ADMno1::Sc()
{
    return Sc_;
}

inline scalar 
Foam::ADMno1::Sct()
{
    return Sct_;
}

inline admPara* 
Foam::ADMno1::paraPtr()
{
    return &para_;
}

inline Foam::PtrList<Foam::volScalarField>& 
Foam::ADMno1::Y()
{
    return YPtrs_;
}

inline Foam::PtrList<Foam::volScalarField::Internal>& 
Foam::ADMno1::dY()
{
    return dYPtrs_;
}

inline PtrList<dimensionedScalar>& 
Foam::ADMno1::GAve()
{
    return GAvePtrs_;
}

inline PtrList<volScalarField>& 
Foam::ADMno1::GR()
{
    return GRPtrs_;
}

inline PtrList<dimensionedScalar>& 
Foam::ADMno1::GRAve()
{
    return GRAvePtrs_;
}

inline PtrList<dimensionedScalar>& 
Foam::ADMno1::dGAve()
{
    return dGAvePtrs_;
}

inline Foam::PtrList<Foam::volScalarField>& 
Foam::ADMno1::M()
{
    return MPtrs_;
}

inline Foam::PtrList<Foam::volScalarField>& 
Foam::ADMno1::IO()
{
    return IOPtrs_;
}

inline Foam::PtrList<Foam::volScalarField::Internal>& 
Foam::ADMno1::dIO()
{
    return dIOPtrs_;
}

// ************************************************************************* //