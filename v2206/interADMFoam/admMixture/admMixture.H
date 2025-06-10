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

// #include "thermoIncompressibleTwoPhaseMixture.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "fvMatrices.H"

#include "fvCFD.H"
#include "upwind.H"
#include "downwind.H"
#include "CMULES.H"
#include "admInterfaceProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class admMixture Declaration
\*---------------------------------------------------------------------------*/

class admMixture
:
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to ADMno1 soluable species
        PtrList<volScalarField>& Si_;

        //- Reference to ADMno1 gaseous species
        PtrList<volScalarField>& Gi_;
        
        // //- Interface properties
        // admInterfaceProperties interface_;

        //- Reference to VoF of liquid phase
        // TODO: maybe not const
        const volScalarField& alpha1_;

        //- VoF of gas phase
        volScalarField alpha2_;

        //- Reference to the thermoIncompressibleTwoPhaseMixture
        // const thermoIncompressibleTwoPhaseMixture& mixture_;

        //- Reference to fvMesh
        const fvMesh& mesh_;

        //- Referemce to flux
        const surfaceScalarField& phi_;

        //- Henry's coefficient
        // TODO: more to come
        dimensionedScalar H_;

        //- Diffusion coeficient of soluable species in liquid phase 
        // TODO: more to come
        dimensionedScalar DS_;

        //- Diffusion coeficient of gaseous species in gas phase
        // TODO: more to come
        dimensionedScalar DG_;

        //- Mean phase diffusion coefficient for soluables in liquid
        surfaceScalarField DalphaS_;

        //- Mean phase diffusion coefficient for gaseous in gas
        surfaceScalarField DalphaG_;

        //- Interface mass transfer flux by Henry's Law
        surfaceScalarField phiHS_;


    // Private Member Functions

        //- No copy construct
        admMixture
        (
            const admMixture&
        ) = delete;

        //- No copy assignment
        void operator=(const admMixture&) = delete;

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
            const admInterfaceProperties& interface
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
                const fvMesh& mesh,
                const volScalarField& alpha1,
                PtrList<volScalarField>& Si,
                PtrList<volScalarField>& Gi,
                const IOdictionary& dict,
                const surfaceScalarField& phi
            ),
            (
                mesh,
                alpha1,
                Si,
                Gi,
                dict,
                phi
            )
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<admMixture> New
        (
            const fvMesh& mesh,
            const volScalarField& alpha1,
            PtrList<volScalarField>& Si,
            PtrList<volScalarField>& Gi,
            const IOdictionary& dict,
            const surfaceScalarField& phi
        );


    // Constructors

        //- Construct from components
        admMixture
        (
            const fvMesh& mesh,
            const volScalarField& alpha1,
            PtrList<volScalarField>& Si,
            PtrList<volScalarField>& Gi,
            const IOdictionary& dict,
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
        void func();

        //- Solve mixture 
        void solve (const admInterfaceProperties& interface);

        inline const volScalarField& alpha() const;

        inline const surfaceScalarField& DalphaS() const;

        inline const surfaceScalarField& DalphaG() const;

        inline const surfaceScalarField& phiHS() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "admMixtureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
