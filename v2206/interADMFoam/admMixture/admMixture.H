/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::admMixture

Description

SourceFiles
    admMixture.C
    newtemperaturePhaseChangeTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef admMixture_H
#define admMixture_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "fvMatrices.H"

#include "fvCFD.H"
#include "upwind.H"
#include "downwind.H"
#include "CMULES.H"
#include "ADMno1.H"
#include "interfaceProperties.H"
#include "incompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class admMixture Declaration
\*---------------------------------------------------------------------------*/

class admMixture
:
    // public IOdictionary
    public incompressibleTwoPhaseMixture
{

protected:

    // Protected data

        // //- Keep a reference to the transportProperties dictionary
        // const dictionary& transportPropertiesDict_;

        //- Surface wettability
        // higher alphaW means smaller contact angle (more wetting)
        dimensionedScalar alphaW_;

        //- Henry's coefficient
        // TODO: maybe remove this eventually?
        dimensionedScalar H_;

        //- Diffusion coeficient of soluable species in liquid phase
        dimensionedScalar D1_;

        //- Diffusion coeficient of gaseous species in gas phase
        dimensionedScalar D2_;

        //- Mean phase diffusion coefficient for soluables in liquid
        surfaceScalarField D1Eff_;

        //- Mean phase diffusion coefficient for gaseous in gas
        surfaceScalarField D2Eff_;

        //- Interface mass transfer flux by Henry's Law
        surfaceScalarField phiHS_;

        // testing ------------------------------------------------------------

        //- ADMno1 reaction
        autoPtr<ADMno1> reaction_;

        //- Phase corrected concentration
        PtrList<volScalarField> SiAlpha_;
        PtrList<volScalarField> GiAlpha_;

        //- Cells marked as mainly fluid
        volScalarField alpha1Full_;
        volScalarField alpha1Empty_;
        volScalarField alpha1Interface_;

        //- Volume fraction threshold to activate inter-phase mass transfer
        dimensionedScalar actAlpha_;

        //- Name of patches to activate inter-phase mass transfer
        wordList actPatch_;

        //- Cells marked as activated for inter-phase mass transfer due to interface
        volScalarField actAlphaCells_;

        //- Cells marked as activated for inter-phase mass transfer due to solid surface
        volScalarField actPatchCells_;

        //- Phase change ADMno1 
        volScalarField mDot_;
        volScalarField mDotAlphal_;
        volScalarField vDot_;
        volScalarField vDotAlphal_;

        // DEBUG scalar
        dimensionedScalar mDotTest_;
        surfaceScalarField phiD_;
        volScalarField Mflux_;
        

    // Private Member Functions

        //- No copy construct
        admMixture
        (
            const admMixture&
        ) = delete;

        //- No copy assignment
        void operator=(const admMixture&) = delete;

        //- Find actAlphaCells_
        void actAlphaCells();
        
        //- Find actPatchCells_
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
        void speciesMules
        (
            const interfaceProperties& interface
        );


public:

    //- Runtime type information
    TypeName("admMixture");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            admMixture,
            components,
            (
                const volScalarField& Top,
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (
                Top,
                U,
                phi
            )
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<admMixture> New
        (
            const volScalarField& Top,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Constructors

        //- Construct from components
        admMixture
        (
            const volScalarField& Top,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~admMixture() = default;


    // Member Functions

        // //- Return the mass condensation and vaporisation rates as a
        // //  coefficient to multiply (1 - alphal) for the condensation rate
        // //  and a coefficient to multiply  alphal for the vaporisation rate
        // virtual Pair<tmp<volScalarField>> mDotAlphal() const = 0;

        // //- Return the mass condensation and vaporisation rates as coefficients
        // virtual Pair<tmp<volScalarField>> mDot() const = 0;

        // //- Return the mass condensation and vaporisation rates as a
        // //  coefficient to multiply (Tsat - T) for the condensation rate
        // //  and a coefficient to multiply  (T - Tsat) for the vaporisation rate
        // virtual Pair<tmp<volScalarField>> mDotDeltaT() const = 0;

        // //- Source for T equarion
        // virtual tmp<fvScalarMatrix> TSource() const = 0;

        // //- Return the volumetric condensation and vaporisation rates as a
        // //  coefficient to multiply (1 - alphal) for the condensation rate
        // //  and a coefficient to multiply  alphal for the vaporisation rate
        // virtual Pair<tmp<volScalarField>> vDotAlphal() const;

        // //- Return the volumetric condensation and vaporisation rates as
        // //  coefficients
        // virtual Pair<tmp<volScalarField>> vDot() const;

        // //- Correct the phaseChange model
        // virtual void correct() = 0;

        // //- Read the transportProperties dictionary and update
        // virtual bool read();

        //- Correct interface
        void limit();

        const volScalarField& mDot();

        const volScalarField& mDotAlphal();

        const volScalarField& vDot();

        const volScalarField& vDotAlphal();

        //- Volume fractions limiter
        void limitAlpha();

        //- Solve mixture 
        void solve
        (
            const interfaceProperties& interface
        );

        //- Solve ADM reactions
        void reaction
        (
            const volScalarField& Top,
            const volScalarField& p
        );

        //- Inline funcs
        inline const volScalarField& alpha() const;

        inline const surfaceScalarField& DalphaS() const;

        inline const surfaceScalarField& DalphaG() const;

        inline const surfaceScalarField& phiHS() const;

        // inline const volScalarField& mDot() const;

        // inline const volScalarField& mDotAlphal() const;

        // inline const volScalarField& vDot() const;

        // inline const volScalarField& vDotAlphal() const;

        inline PtrList<volScalarField>& Si();

        inline PtrList<volScalarField>& Gi();

        inline PtrList<volScalarField>& SiAlpha();

        inline PtrList<volScalarField>& GiAlpha();

        // testing ------------------------------------------------------------
        // TODO: corrected by mesh().V()
        dimensionedScalar speciesAlphaAverage
        (
            const volScalarField& YiAlpha,
            const volScalarField& alpha
        )
        {
            volScalarField magAlpha = mag(alpha);

            return
            (
                dimensionedScalar
                (
                    YiAlpha.dimensions()/alpha.dimensions(),
                    gSum(YiAlpha)/gSum(magAlpha)
                )
            );
        }


        void speciesAlphaCorrect()
        {
            forAll(SiAlpha_, i)
            {
                SiAlpha_[i] = reaction_->Y()[i]*alpha1_;
            }

            forAll(GiAlpha_, i)
            {
                GiAlpha_[i] = reaction_->G()[i]*alpha2_;
            }

            // TODO: for source terms dS and dG too
        }


        void speciesADMCorrect()
        {
            //- threshold volume fraction alpha0, to avoid error introduced by interface correction
            // alpha0 = 0.10  -->  high error
            // alpha0 = 0.90  -->  low error
            // limitAlpha();
            scalar alpha0 = 0.9;
            scalar alpha00 = 1e-2;

            dimensionedScalar YiAlphaAverage
            (
                dimMass/dimVolume,
                Zero
            );

            dimensionedScalar YiAlphaAverageInterface
            (
                dimMass/dimVolume,
                Zero
            );

            volScalarField limitedAlpha1
            (
                min(max(alpha1_, scalar(1e-16)), scalar(1))
            );

            // forAll(SiAlpha_, i)
            // {
                label i = 6;

                alpha1Empty_ = min(zeroField(), limitedAlpha1 - alpha00)/(limitedAlpha1 - alpha00);
                alpha1Full_ = max(zeroField(), limitedAlpha1 - alpha0)/(limitedAlpha1 - alpha0);
                alpha1Interface_ = 1 - alpha1Empty_ - alpha1Full_;

                // Total mass and volume field
                volScalarField::Internal massField = SiAlpha_[i].internalField()*U_.mesh().V();
                volScalarField::Internal volumeField = limitedAlpha1*alpha1Full_.mesh().V();
                YiAlphaAverage.value() = gSum(massField) / gSum(volumeField);

                // // Total mass from full and full cells
                // volScalarField::Internal massFieldFull = alpha1Full_*SiAlpha_[i].internalField()*alpha1Full_*U_.mesh().V();
                // volScalarField::Internal volumeFieldFull = alpha1Full_*limitedAlpha1*alpha1Full_.mesh().V();

                // Mass and volume field at interface cells
                volScalarField::Internal massFieldInterface = alpha1Interface_*SiAlpha_[i].internalField()*U_.mesh().V();
                volScalarField::Internal volumeFieldInterface = alpha1Interface_*limitedAlpha1.internalField()*U_.mesh().V();
                if (gSum(alpha1Interface_) < 1)
                {
                    YiAlphaAverageInterface.value() = YiAlphaAverage.value();
                }
                else
                {
                    YiAlphaAverageInterface.value() = gSum(massFieldInterface) / gSum(volumeFieldInterface);   
                }

                // Recreate ADM Yi fields
                reaction_->Y()[i] = 
                (
                    SiAlpha_[i]/limitedAlpha1*alpha1Full_ 
                  + YiAlphaAverage*alpha1Empty_
                  + YiAlphaAverageInterface*alpha1Interface_
                );

                // DEBUG
                volScalarField::Internal massFieldADMAlpha = reaction_->Y()[i].internalField()*limitedAlpha1*alpha1Full_.mesh().V();
                volScalarField::Internal massFieldADMcorrected = reaction_->Y()[i].internalField()*0.75*alpha1Full_.mesh().V();
                volScalarField::Internal volumeFieldADM = limitedAlpha1*alpha1Full_.mesh().V();

                Info<< "DEBUG!!!" << endl;
                volScalarField::Internal massFieldADM = reaction_->Y()[i].internalField();//*limitedAlpha1*alpha1Full_.mesh().V();
                Info<<   ">>> Sac mass (* alpha): " << gSum(massFieldADMAlpha)
                    << "\n>>> Sac mass (w/o alpha): " << gSum(massFieldADMcorrected)
                    << "\n>>> phase volume: " << gSum(volumeFieldADM) << endl; 
                // Info<< ">>> " << gSum(alphaAll_) 
                //     << " = " << gSum(alpha1Empty_) + gSum(alpha1Full_)
                //     << " ? " << gSum(alpha1Empty_) + gSum(alpha1Full_) + gSum(alpha1Interface_) << endl;

            // // -------------------------------------------------------------------------------------------------
            //     volScalarField limitedAlpha1
            //     (
            //         min(max(alpha1_, scalar(1e-6)), scalar(1))
            //     );

            //     alpha1Full_ = max(zeroField(), limitedAlpha1 - alpha0)/(limitedAlpha1 - alpha0);
            //     volScalarField alpha1Empty = 1 - alpha1Full_;

            //     // total mass of species in non-full cells
            //     volScalarField::Internal YiMassEmpty = SiAlpha_[i].internalField()*alpha1Empty*alpha1_.mesh().V();
            //     dimensionedScalar residualYiMassEmpty = gSum(YiMassEmpty);

            //     // total volume of non-full cells
            //     volScalarField::Internal alpha1EmptyV = alpha1Empty*alpha1_.mesh().V();
            //     dimensionedScalar residualV = gSum(alpha1EmptyV);

            //     // average volume fraction of alpha1 in non-full cells
            //     volScalarField::Internal alpha1Alpha1Empty = alpha1Empty*limitedAlpha1;
            //     dimensionedScalar residualAlpha = gSum(alpha1Alpha1Empty)/gSum(alpha1Empty);

            //     //
            //     YiAlphaAverage.value() = (residualYiMassEmpty/residualAlpha/residualV).value();
            //     reaction_->Y()[i] = SiAlpha_[i]/limitedAlpha1*alpha1Full_ + YiAlphaAverage*alpha1Empty;

            //     //
            //     volScalarField::Internal sum = reaction_->Y()[i].internalField()*limitedAlpha1.internalField()*alpha1_.mesh().V();
            //     Info<< "Species concentration (ADM corrected sum) = "
            //         << gSum(sum) << endl;
            //     Info<< ">>> total mass of species in non-full cells = "
            //         << residualYiMassEmpty.value() << endl;
            //     Info<< ">>> total volume of non-full cells = "
            //         << residualV.value() << endl;
            //     Info<< ">>> average volume fraction of alpha1 in non-full cells = "
            //         << residualAlpha.value() << endl;
            //     Info<< ">>> YiAlphaAverage for empty = "
            //         << YiAlphaAverage.value() << endl;
            // // -------------------------------------------------------------------------------------------------

            //     // YiAlphaAverage = speciesAlphaAverage(SiAlpha_[i], alpha1_);
            //     // reaction_->Y()[i] = SiAlpha_[i]/alpha1_*alpha1Full_ + YiAlphaAverage*(1 - alpha1Full_);
            // // }

            // forAll(GiAlpha_, i)
            // {
            //     alpha1Full_ = max(zeroField(), alpha2_ - alpha0)/(alpha2_ - alpha0);
            //     YiAlphaAverage = speciesAlphaAverage(GiAlpha_[i], alpha2_);

            //     reaction_->G()[i] = GiAlpha_[i]/alpha2_*alpha1Full_ + YiAlphaAverage*(1 - alpha1Full_);
            // }
        }


        void printGasGenRate()
        {
            // total gas generation rate
            // dimensionedScalar sumGasGen = - gSum(mDot_);
            // sumGasGen.dimensions().reset(mDot_.dimensions());
            Info<< "total gas generation rate [mol * m-3]: " << - mDot_.weightedAverage(alpha1_.mesh().V()).value() << endl;

            // forAll(mDot_, i)
            // {
            //     dimensionedScalar sumGasGen = gSum(mDot_[i]);
            //     Info<< "total gas generation rate [mol * m-3]: " << sumGasGen << endl;
            // }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "admMixtureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
