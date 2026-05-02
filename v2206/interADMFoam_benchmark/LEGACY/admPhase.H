/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::multiphaseADMixture

Description
    Incompressible multi-phaseADM mixture with built in solution for the
    phaseADM fractions with interface compression for interface-capturing.

    Derived from transportModel so that it can be unused in conjunction with
    the incompressible turbulence models.

    Surface tension and contact-angle is handled for the interface
    between each phaseADM-pair.

SourceFiles
    multiphaseADMixture.C

\*---------------------------------------------------------------------------*/

#ifndef admPhase_H
#define admPhase_H

// #include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"
#include "incompressible/transportModel/transportModel.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "PtrDictionary.H"
#include "PtrListDictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrix.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations


/*---------------------------------------------------------------------------*\
                      Class multiphaseADMixture Declaration
\*---------------------------------------------------------------------------*/

class admPhase
:
    public IOdictionary,
    public transportModel
{
protected:
 
    // Protected data

        word phase1Name_;
        word phase2Name_;
 
        volScalarField alpha1_;
        volScalarField alpha2_;
 
        autoPtr<viscosityModel> nuModel1_;
        autoPtr<viscosityModel> nuModel2_;
 
        dimensionedScalar rho1_;
        dimensionedScalar rho2_;
 
        const volVectorField& U_;
        const surfaceScalarField& phi_;
 
        volScalarField nu_;
 
 
    // Protected Member Functions
 
        //- Calculate and return the laminar viscosity
        void calcNu();
 
 
public:
 
    TypeName("admPhase");
 
 
    // Constructors
 
        //- Construct from components
        admPhase
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );
 
 
    //- Destructor
    virtual ~admPhase() = default;
 
 
    // Member Functions

        const word& phase1Name() const
        {
            return phase1Name_;
        }
 
        const word& phase2Name() const
        {
            return phase2Name_;
        }
 
        //- Return the phase-fraction of phase 1
        const volScalarField& alpha1() const
        {
            return alpha1_;
        }
 
        //- Return the phase-fraction of phase 1
        volScalarField& alpha1()
        {
            return alpha1_;
        }
 
        //- Return the phase-fraction of phase 2
        const volScalarField& alpha2() const
        {
            return alpha2_;
        }
 
        //- Return the phase-fraction of phase 2
        volScalarField& alpha2()
        {
            return alpha2_;
        }

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return *nuModel1_;
        }
 
        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return *nuModel2_;
        }
 
        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }
 
        //- Return const-access to phase2 density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };
 
        //- Return const-access to the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }
 
        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;
 
        //- Return the dynamic laminar viscosity on patch
        tmp<scalarField> mu(const label patchI) const;
 
        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;
 
        //- Return the kinematic laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }
 
        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }
 
        //- Return the face-interpolated kinematic laminar viscosity
        tmp<surfaceScalarField> nuf() const;
 
        //- Correct the laminar viscosity
        virtual void correct()
        {
            calcNu();
        }
 
        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //