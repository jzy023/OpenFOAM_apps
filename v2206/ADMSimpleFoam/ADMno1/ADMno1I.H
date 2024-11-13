/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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

inline scalar Foam::ADMno1::Sc()
{
    return Sc_;
}

inline admPara* Foam::ADMno1::paraPtr()
{
    return &para_;
}

inline Foam::PtrList<Foam::volScalarField>& Foam::ADMno1::Y()
{
    return YPtrs_;
}

inline Foam::PtrList<Foam::volScalarField::Internal>& Foam::ADMno1::dY()
{
    return dYPtrs_;
}

inline Foam::PtrList<Foam::volScalarField>& Foam::ADMno1::G()
{
    return GPtrs_;
}

inline Foam::PtrList<Foam::volScalarField::Internal>& Foam::ADMno1::dG()
{
    return dGPtrs_;
}

inline Foam::PtrList<Foam::volScalarField>& Foam::ADMno1::IO()
{
    return IOPtrs_;
}

inline Foam::PtrList<Foam::volScalarField::Internal>& Foam::ADMno1::dIO()
{
    return dIOPtrs_;
}

// inline PtrList<volScalarField>& ADMno1::Y()
// {
//     return YPtrs_;
// }

// inline const Foam::PtrList<Foam::volScalarField>& Y() const
// {
//     return YPtrs_;
// }

// ************************************************************************* //