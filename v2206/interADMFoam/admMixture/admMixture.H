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

        //- ADMno1 reaction
        autoPtr<ADMno1> reaction_;

        // //- Keep a reference to the transportProperties dictionary
        // const dictionary& transportPropertiesDict_;

        //- Surface wettability
        // higher alphaW means smaller contact angle (more wetting)
        dimensionedScalar alphaW_;

        //- Ideal gas constant
        dimensionedScalar R_;

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

        //- Phase corrected concentration
        PtrList<volScalarField> SiAlpha_;
        PtrList<volScalarField> GiAlpha_;

        //- Volume fraction threshold to activate inter-phase mass transfer
        dimensionedScalar alphaInterface_;

        //- Name of patches to activate inter-phase mass transfer
        wordList actPatch_;

        //- Cells marked as mainly fluid
        // TODO: use Internal? and .ref()
        volScalarField isCellsFull_;
        volScalarField isCellsEmpty_;

        //- Cells marked as activated for inter-phase mass transfer due to interface
        volScalarField isCellsInterface_;

        //- Cells marked as activated for inter-phase mass transfer due to solid surface
        volScalarField isCellsActWall_;

        //- Cells marked as activated for overall inter-phase mass transfer 
        volScalarField isCellsAct_;

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

        //- Find isCellsAct_
        void findCellAct();
        
        //- Find isCellsActWall_
        void findCellsActWall();

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
        void solvePhase
        (
            const interfaceProperties& interface
        );

        //- Solve ADM reactions
        void solveReaction
        (
            const surfaceScalarField &flux,
            const volScalarField& Top
            // const volScalarField& Top,
            // const volScalarField& p
        );

        //- Update Gas Concentration
        void updateReactionGas
        (
            const dimensionedScalar deltaT
        );

        //- Inline functions
        inline const volScalarField& alpha() const;

        inline const volScalarField& actCells() const;

        inline const surfaceScalarField& DalphaS() const;

        inline const surfaceScalarField& DalphaG() const;

        inline const surfaceScalarField& phiHS() const;

        // inline const volScalarField& mDot() const;

        // inline const volScalarField& mDotAlphal() const;

        // inline const volScalarField& vDot() const;

        // inline const volScalarField& vDotAlphal() const;

        inline PtrList<volScalarField::Internal>& dSi();

        inline PtrList<volScalarField::Internal>& dGi();

        inline PtrList<volScalarField>& SiAlpha();

        inline PtrList<volScalarField>& GiAlpha();

        // testing ------------------------------------------------------------
        // TODO: corrected by mesh().V()
        // dimensionedScalar speciesAlphaAverage
        // (
        //     const volScalarField& YiAlpha,
        //     const volScalarField& alpha
        // )
        // {
        //     volScalarField magAlpha = mag(alpha);

        //     return
        //     (
        //         dimensionedScalar
        //         (
        //             YiAlpha.dimensions()/alpha.dimensions(),
        //             gSum(YiAlpha)/gSum(magAlpha)
        //         )
        //     );
        // }


        void gasPartialPressure
        (
            const volScalarField& Top
        )
        {
            
            // Identifying region for averaging gaseous species for partial pressure calculation
            volScalarField limitedAlpha2
            (
                min(max((1 - alpha1_), scalar(1e-16)), scalar(1))
            );

            volScalarField alphaGasRegion_ = max(zeroField(), limitedAlpha2 - 0.1)/(limitedAlpha2 - 0.1);
            volScalarField::Internal GiAlphaTopKH = 
            (
                alphaGasRegion_.internalField()*GiAlpha_[0].internalField()
              * Top.internalField()*reaction_->KHh2()
            );

            dimensionedScalar GiAlphaTopKHAve = GiAlphaTopKH.weightedAverage(alpha1_.mesh().V());
            dimensionedScalar Ph2 = R_ * GiAlphaTopKHAve;

            // // DEBUG
            // Info<< ">>> " << Ph2.dimensions() << endl;
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

            forAll(SiAlpha_, i)
            {
                volScalarField& Yi = SiAlpha_[i];
                volScalarField& YiADM = reaction_->Y()[i];

                isCellsEmpty_ = min(zeroField(), limitedAlpha1 - alpha00)/(limitedAlpha1 - alpha00);
                isCellsFull_ = max(zeroField(), limitedAlpha1 - (1.0 - alphaInterface_))/(limitedAlpha1 - (1.0 - alphaInterface_));
                isCellsInterface_ = 1 - isCellsEmpty_ - isCellsFull_;

                // Total mass and volume field
                volScalarField::Internal massField = Yi.internalField()*U_.mesh().V();
                volScalarField::Internal volumeField = limitedAlpha1*isCellsFull_.mesh().V();
                YiAlphaAverage.value() = gSum(massField) / gSum(volumeField);

                // // Total mass from full and full cells
                // volScalarField::Internal massFieldFull = isCellsFull_*Yi.internalField()*isCellsFull_*U_.mesh().V();
                // volScalarField::Internal volumeFieldFull = isCellsFull_*limitedAlpha1*isCellsFull_.mesh().V();

                // Mass and volume field at interface cells
                volScalarField::Internal massFieldInterface = isCellsInterface_*Yi.internalField()*U_.mesh().V();
                volScalarField::Internal volumeFieldInterface = isCellsInterface_*limitedAlpha1.internalField()*U_.mesh().V();
                if (gSum(isCellsInterface_) < 1)
                {
                    YiAlphaAverageInterface.value() = YiAlphaAverage.value();
                }
                else
                {
                    YiAlphaAverageInterface.value() = gSum(massFieldInterface) / gSum(volumeFieldInterface);   
                }

                // Recreate ADM Yi fields
                YiADM = 
                (
                    Yi/limitedAlpha1*isCellsFull_ 
                  + YiAlphaAverage*isCellsEmpty_
                  + YiAlphaAverageInterface*isCellsInterface_
                );

                // DEBUG
                volScalarField::Internal massFieldADMAlpha = YiADM.internalField()*limitedAlpha1*isCellsFull_.mesh().V();
                volScalarField::Internal massFieldADMcorrected = YiADM.internalField()*0.75*isCellsFull_.mesh().V();
                volScalarField::Internal volumeFieldADM = limitedAlpha1*isCellsFull_.mesh().V();

                // Info<< "DEBUG!!!" << endl;
                // volScalarField::Internal massFieldADM = YiADM.internalField();//*limitedAlpha1*isCellsFull_.mesh().V();
                // Info<<   ">>> Sac mass (* alpha): " << gSum(massFieldADMAlpha)
                //     << "\n>>> Sac mass (w/o alpha): " << gSum(massFieldADMcorrected)
                //     << "\n>>> phase volume: " << gSum(volumeFieldADM) << endl;
            }
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

        // autoPtr<ADMno1> reaction()
        ADMno1* reaction()
        {
            return reaction_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "admMixtureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
