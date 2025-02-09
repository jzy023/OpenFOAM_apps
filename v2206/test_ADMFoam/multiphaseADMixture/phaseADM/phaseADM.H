/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::phaseADM

Description
    Single incompressible phaseADM derived from the phaseADM-fraction.
    Used as part of the multiphaseADMixture for interface-capturing multi-phaseADM
    simulations.

SourceFiles
    phaseADM.C

\*---------------------------------------------------------------------------*/

#ifndef phaseADM_H
#define phaseADM_H

#include "volFields.H"
#include "dictionaryEntry.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseADM Declaration
\*---------------------------------------------------------------------------*/

class phaseADM
:
    public volScalarField
{
    // Private data

        word name_;
        dictionary phaseDict_;
        autoPtr<viscosityModel> nuModel_;
        
        // testing
        bool isRhoField_ = false;
        dimensionedScalar rho_;


public:

    // Constructors

        //- Construct from components
        phaseADM
        (
            const word& name,
            const dictionary& phaseDict,
            const volVectorField& U,
            const surfaceScalarField& phi
        );

        //- Return clone
        autoPtr<phaseADM> clone() const;

        //- Return a pointer to a new phaseADM created on freestore
        //  from Istream
        class iNew
        {
            const volVectorField& U_;
            const surfaceScalarField& phi_;

        public:

            iNew
            (
                const volVectorField& U,
                const surfaceScalarField& phi
            )
            :
                U_(U),
                phi_(phi)
            {}

            autoPtr<phaseADM> operator()(Istream& is) const
            {
                dictionaryEntry ent(dictionary::null, is);
                return autoPtr<phaseADM>(new phaseADM(ent.keyword(), ent, U_, phi_));
            }
        };


    // Member Functions

        const word& name() const
        {
            return name_;
        }

        const word& keyword() const
        {
            return name();
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

        //- Return const-access to phase1 density
        const dimensionedScalar& rho() const
        {
            return rho_;
        }

        //- Correct the phaseADM properties
        void correct();

        //-Inherit read from volScalarField
        using volScalarField::read;

        //- Read base transportProperties dictionary
        bool read(const dictionary& phaseDict);


    // testing
        //- Set rho
        scalar setRho()
        {
            if (isRhoField_)
            {
                Info<< "Setting density field to be calculated" 
                    << " for phase: " << name_ << endl;
                return 1.0; // dummy value, rho filed would be calculated in mixture
            }
            else
            {
                Info<< "Setting density field to be homogeneous" 
                    << " for phase: " << name_ << endl;
                return phaseDict_.getScalar("rho");
            }
        }

        //- Set rho
        bool setIsRhoField()
        {
            if (name_ == "liquid")
            {
                return false;
            }
            else
            {
                bool ans =
                (
                    !phaseDict_.found("rho") | 
                    phaseDict_.lookupOrDefault("isRhoField", false)
                );
                
                return ans;
            }
        }

        //- return isRhoField_  
        bool isRhoField() const
        {  
            return isRhoField_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
