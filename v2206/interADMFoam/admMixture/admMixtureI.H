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

inline const Foam::volScalarField&
Foam::admMixture::alpha() const
{
    return alpha1_;
}


inline const Foam::volScalarField&
Foam::admMixture::actCells() const
{
    return isCellsAct_;
}


inline const Foam::surfaceScalarField&
Foam::admMixture::DalphaS() const
{
    return D1Eff_;
}


inline const Foam::surfaceScalarField&
Foam::admMixture::DalphaG() const
{
    return D2Eff_;
}


inline const Foam::surfaceScalarField&
Foam::admMixture::phiHS() const
{
    return phiHS_;
}


// inline const Foam::volScalarField&
// Foam::admMixture::mDot() const
// {
//     return mDot_;
// }


// inline const Foam::volScalarField&
// Foam::admMixture::mDotAlphal() const
// {
//     return mDotAlphal_;
// }


// inline const Foam::volScalarField&
// Foam::admMixture::vDot() const
// {
//     return vDot_;
// }


// inline const Foam::volScalarField&
// Foam::admMixture::vDotAlphal() const
// {
//     return vDotAlphal_;
// }


inline Foam::PtrList<Foam::volScalarField::Internal>&
Foam::admMixture::dSi()
{
    return reaction_->dY();
}


inline Foam::PtrList<Foam::volScalarField::Internal>&
Foam::admMixture::dGi()
{
    return reaction_->dG();
}


inline Foam::PtrList<Foam::volScalarField>&
Foam::admMixture::SiAlpha()
{
    return SiAlpha_;
}


inline Foam::PtrList<Foam::volScalarField>&
Foam::admMixture::GiAlpha()
{
    return GiAlpha_;
}


// ************************************************************************* //