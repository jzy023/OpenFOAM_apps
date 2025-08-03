/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::admPhases
 
Description
    A two-phase model for ADMno1
 
SourceFiles
    admPhases.C
 
\*---------------------------------------------------------------------------*/
 
#ifndef admPhases_H
#define admPhases_H
 
#include "volFields.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
                      Class admPhases Declaration
\*---------------------------------------------------------------------------*/
 
class admPhases
:
    public volScalarField
{
protected:
 
    // Protected data
 
        // volScalarField alpha_;
        
        //- Phase name
        word name_;

        //- Phase density
        dimensionedScalar rho_;

        //- Phase visocsity model
        autoPtr<viscosityModel> nuModel_;

        //- Phase volumetric flux
        surfaceScalarField alphaPhi_;
 
 
public:
 
    // Constructors
 
        //- Construct from components
        admPhases
        (
            const fvMesh& mesh,
            const dictionary& dict
        );
 
 
    //- Destructor
    ~admPhases()
    {}
 
 
    // Member Functions

        //- Return the name of phase 1
        const word& name() const
        {
            return name_;
        }

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel() const
        {
            return *nuModel_;
        }

        //- Return the kinematic laminar viscosity
        tmp<volScalarField> nu() const
        {
            return nuModel_->nu();
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return nuModel_->nu(patchi);
        }

        //- Return phase density
        dimensionedScalar& rho() const
        {
            return rho_;
        }

        //- Return const-access to phase density
        const dimensionedScalar& rho() const
        {
            return rho_;
        }

        //- Return phase volumetric flux
        surfaceScalarField& alphaPhi()
        {
            return alphaPhi_;
        }

        //- Return const-access to phase volumetric flux
        const surfaceScalarField& alphaPhi() const
        {
            return alphaPhi_;
        }
};
 
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
} // End namespace Foam
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#endif
 
// ************************************************************************* //