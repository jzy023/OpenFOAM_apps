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

#include "ADMno1Types.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    inline scalar admPara::DTOS() const 
    {
        return DTOS_;
    }

    inline scalar admPara::MTOm() const 
    {
        return MTOm_;
    }

    inline scalar admPara::BTOP() const 
    {
        return BTOP_;
    }

    inline scalar admPara::kTOK() const 
    {
        return kTOK_;
    }

    inline const dimensionedScalar admPara::Tbase() const 
    {
        return Tbase_;
    }

    inline const dimensionedScalar admPara::T() const 
    {
        return Top_;
    }

    inline const decayRate admPara::kDec() const 
    {
        return kDec_;
    }

    inline const yieldBiomass admPara::yB() const 
    {
        return yB_;
    }

    inline const yieldProduct admPara::yP() const 
    {
        return yP_;
    }

    inline const carbonContent admPara::CC() const 
    {
        return CC_;
    }

    inline const inhibitionParaI admPara::KI() const 
    {
        return KI_;
    }

    inline const inhibitionParaS admPara::KS() const 
    {
        return KS_;
    }

    inline const admAcidBaseEq admPara::Ka() const 
    {
        return Ka_;
    }

    inline const gasHenryCoe admPara::KH() const 
    {
        return KH_;
    }

    inline const dimensionedScalar admPara::kLa() const 
    {
        return kLa_;
    }

    inline const limPH admPara::pHL() const 
    {
        return pHL_;
    }

    inline const nitrogenContent admPara::NC() const 
    {
        return NC_;
    }

    inline scalar admPara::Gini
    (
        label idx
    ) const
    {
        return Gini_[idx]; 
    }

    inline scalar admPara::Mini
    (
        label idx
    ) const
    {
        return Mini_[idx]; 
    }

    inline scalar admPara::Pini() const 
    {
        return Pini_;
    }

    inline scalar admPara::Eini
    (
        label idx
    ) const
    {
        return Eini_[idx]; 
    }

    inline scalar admPara::Iini
    (
        label idx
    ) const
    {
        return Iini_[idx]; 
    }

// ************************************************************************* //