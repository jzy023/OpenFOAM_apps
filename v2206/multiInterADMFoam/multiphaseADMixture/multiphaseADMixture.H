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

#ifndef multiphaseADMixture_H
#define multiphaseADMixture_H

#include "incompressible/transportModel/transportModel.H"
#include "IOdictionary.H"
#include "PtrDictionary.H"
#include "PtrListDictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrix.H"

#include "phaseADM.H"
#include "ADMno1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class alphaContactAngleFvPatchScalarField;

/*---------------------------------------------------------------------------*\
                      Class multiphaseADMixture Declaration
\*---------------------------------------------------------------------------*/

class multiphaseADMixture
:
    public IOdictionary,
    public transportModel
{
public:

    //- Symmetric pair of interface names
    class interfacePair
    :
        public Pair<word>
    {
    public:

        // Always use symmetric hashing
        using hasher = Pair<word>::symmHasher;

        // Always use symmetric hashing (alias)
        using hash = Pair<word>::symmHasher;


        // Constructors

            interfacePair() = default;

            interfacePair(const word& alpha1Name, const word& alpha2Name)
            :
                Pair<word>(alpha1Name, alpha2Name)
            {}

            interfacePair(const phaseADM& alpha1, const phaseADM& alpha2)
            :
                Pair<word>(alpha1.name(), alpha2.name())
            {}


        // Friend Operators

            friend bool operator==
            (
                const interfacePair& a,
                const interfacePair& b
            )
            {
                return (0 != Pair<word>::compare(a, b));
            }

            friend bool operator!=
            (
                const interfacePair& a,
                const interfacePair& b
            )
            {
                return (!(a == b));
            }
    };


private:

    // Private data

        //- Dictionary of phases
        PtrDictionary<phaseADM> phases_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to the velocity field
        const volVectorField& U_;

        //- Reference to the volumetric flux
        const surfaceScalarField& phi_;

        //- Mixture total mass flux
        surfaceScalarField rhoPhi_;

        //- Volume fractions of each phase
        volScalarField alphas_;

        //- Kinetic viscocity
        volScalarField nu_;

        //- Data type for table of surface tension coefficients
        typedef HashTable<scalar, interfacePair, interfacePair::hash>
            sigmaTable;

        //- Table of surface tension coefficients
        sigmaTable sigmas_;

        //- Surface tension coefficients dimensions
        dimensionSet dimSigma_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        //- Interface normal surface vector
        surfaceScalarField nHatf_;

        //- ADMno1 reaction
        autoPtr<ADMno1> reaction_;

        //- Mass transfer phase source terms
        typedef HashTable<volScalarField::Internal> SuSpTable;

        //- Su phase source terms
        SuSpTable Su_;
 
        //- Sp phase source terms
        SuSpTable Sp_;

        //- Maximum volumen rate change
        dimensionedScalar ddtAlphaMax_;

        // TODO: future implementations
        // //- Turbulent Prandt number
        // dimensionedScalar Prt_;

        // testing -----------------------------------------------------------
        //- Phase corrected concentration
        PtrList<volScalarField> SiAlpha_;
        PtrList<volScalarField> GiAlpha_;

        //- Cells marked as mainly fluid
        // TODO: rename
        volScalarField alpha1Full_;

        //- Volume fraction threshold to activate inter-phase mass transfer
        dimensionedScalar actAlpha_;

        //- Name of patches to activate inter-phase mass transfer
        wordList actPatch_;

        //- Cells marked as activated for inter-phase mass transfer due to interface
        volScalarField::Internal actGasMixCells_;

        //- Cells marked as activated for inter-phase mass transfer due to solid surface
        volScalarField::Internal actPatchCells_;

        // //- Phase change ADMno1 
        // volScalarField mDot_;
        // volScalarField mDotAlphal_;
        // volScalarField vDot_;
        // volScalarField vDotAlphal_;

        //- Surface wettability
        // higher alphaW means smaller contact angle (more wetting)
        dimensionedScalar alphaW_;
        
        // DEBUG scalar
        dimensionedScalar mDotTest_;
        // surfaceScalarField phiD_;
        // volScalarField Mflux_;
        // -------------------------------------------------------------------


    // Private member functions

        //- 
        void calcAlphas();

        //- Solve alphas
        void solveAlphas
        (
            const scalar cAlpha
            // ,const PtrListDictionary<volScalarField::Internal>& vDotList
        );

        //- Interface normal surface vector
        tmp<surfaceVectorField> nHatfv
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Calculate and update interface normal surface vector
        tmp<surfaceScalarField> nHatf
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        );

        //- Correction for nHat 
        void correctContactAngle
        (
            const phaseADM& alpha1,
            const phaseADM& alpha2,
            surfaceVectorField::Boundary& nHatb
        ) const;

        //- Correction for the boundary condition for nHat on walls
        void correctBoundaryContactAngle
        (
            const alphaContactAngleFvPatchScalarField& acap,
            label patchi,
            const phaseADM& alpha1,
            const phaseADM& alpha2,
            surfaceVectorField::Boundary& nHatb
        ) const;

        //- Interface Curvature
        tmp<volScalarField> K(const phaseADM& alpha1, const phaseADM& alpha2) const;

        // testing -----------------------------------------------------------
        //- Find liquid-gas interface cells
        void actGasMixCells();
        
        //- Find activated wall cells 
        void actPatchCells();

        //- Normalizing for MULES
        void normalizeYi();

        //- Species transport equations
        void massTransferCoeffs();

        //- Interface compression coefficient 
        surfaceScalarField compressionCoeff
        (
            const volScalarField& Yi
        );

        //- MULES solver for species
        void speciesMules();
        // -------------------------------------------------------------------


public:

    // Constructors

        //- Construct from components
        multiphaseADMixture
        (
            const volScalarField& Top,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~multiphaseADMixture() = default;


    // Member Functions

        //- Return const-access to the phases
        const PtrDictionary<phaseADM>& phases() const
        {
            return phases_;
        }

        //- Return the velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return const-access to the volumetric flux
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- Return const-access to the mass flux
        const surfaceScalarField& rhoPhi() const
        {
            return rhoPhi_;
        }

        //- Access to ddtAlphaMax
        dimensionedScalar ddtAlphaMax() const
        {
            return ddtAlphaMax_;
        }

        //- Return the mixture density
        tmp<volScalarField> rho() const;

        //- Return the mixture density for patch
        tmp<scalarField> rho(const label patchi) const;

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;

        //- Return the dynamic laminar viscosity for patch
        tmp<scalarField> mu(const label patchi) const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- Return the kinematic laminar viscosity
        tmp<volScalarField> nu() const;

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        tmp<surfaceScalarField> surfaceTensionForce() const;

        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;

        //- Source term for pEqn due to ADMno1 gas generation
        tmp<fvScalarMatrix> volTransfer
        (
            const volScalarField& p
            // ,const PtrListDictionary<volScalarField::Internal>& vDotList
        );

        //- Solve for the mixture phaseADM-fractions
        void solve
        (
            // const PtrListDictionary<volScalarField::Internal>& vDotList
        );

        //- Correct the mixture properties
        void correct();

        //- Read base transportProperties dictionary
        bool read();


    // testing -----------------------------------------------------------

        //- Safe check for setting rho
        void checkPhases();

        //- Rerturn name of phases
        const PtrList<word> phaseNames()
        {
            PtrList<word> phaseNames(phases_.size());
            
            label phasei = 0;
            for (phaseADM& phase : phases_)
            {
                phaseNames.set
                (
                    phasei,
                    new word
                    (
                        phase.name()
                    )
                );

                ++phasei;
            }

            return phaseNames;
        }

        //- Return const-access to the interface normal surface vector
        const surfaceScalarField& nHatf() const 
        {
            return nHatf_;
        }

        //- Maximum diffusion number
        scalar maxDiffNo() const
        {   // TODO: fix?
            return 0.0;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
