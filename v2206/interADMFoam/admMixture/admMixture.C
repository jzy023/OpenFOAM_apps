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

\*---------------------------------------------------------------------------*/

#include "admMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admMixture, 0);
    defineRunTimeSelectionTable
    (
        admMixture,
        components
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Species transport equations
void Foam::admMixture::massTransferCoeffs()
{
    alpha2_ = 1.0 - alpha1_;
            
    // calculate and return mean diffusion coefficient
    // TODO: multispecies? add turbulent diffusivity too?
    DalphaG_ = fvc::interpolate(DG_ * alpha2_);
    DalphaS_ = 
    (
        fvc::interpolate(DS_ * alpha1_)
        // fvc::interpolate(DS_ * alpha1_ + H_ * DS2_ * alpha2_) / fvc::interpolate(alpha1_ + H_ * alpha2_)
    );
            
    //- Calculate interface mass transfer flux by Henry's Law
    // TODO: check for gas [all 0 if H_ -> 0]
    // surfaceScalarField phiHUp = speciesMixture.phiHUp(i);
    // surfaceScalarField phiHDown = speciesMixture.phiHDown(i);
    phiHS_ = 
    (
        DalphaS_ * (1 - H_) / fvc::interpolate((alpha1_ + H_ * (1 - alpha1_)))
      * fvc::snGrad(alpha1_) * U_.mesh().magSf()
    );
}


//- interface compreassion coefficient
surfaceScalarField Foam::admMixture::compressionCoeff
(
    const volScalarField& Yi
)
{
    // Direction of interfacial flux
    surfaceScalarField fluxDir = fvc::snGrad(alpha1_) * U_.mesh().magSf();

    // Upwind and downwind alpha1
    surfaceScalarField alphaUp = upwind<scalar>(U_.mesh(),fluxDir).interpolate(alpha1_);
    surfaceScalarField alphaDown = downwind<scalar>(U_.mesh(),fluxDir).interpolate(alpha1_);

    // Upwind and downwnd Yi
    surfaceScalarField YiUp = upwind<scalar>(U_.mesh(),fluxDir).interpolate(Yi);
    surfaceScalarField YiDown = downwind<scalar>(U_.mesh(),fluxDir).interpolate(Yi);
        
    dimensionedScalar sgn = 
    (
        sign(max(alphaDown * YiDown) - max((1 - alphaUp) * YiUp))
    );

    // Normal compression coefficient
    surfaceScalarField deltaYi1 = 
    (
        max
        (
          - max(Yi),
            min
            (
                max(Yi),
                (YiDown - YiUp) / (alphaDown - alphaUp + 1e-4)
            )
        )
    );
        
    // Standard compression coefficient
    surfaceScalarField deltaYi2 = 
    (
        max
        (
          - max(Yi),
            min
            (
                max(Yi),
                (
                    YiDown / (alphaDown + (1 - alphaDown) * H_)
                  - H_ * YiUp / (alphaUp + (1 - alphaUp) * H_)
                )
            )
        )
    );

    return sgn * max(mag(deltaYi1),mag(deltaYi2));
}


void Foam::admMixture::speciesMules
(
    const interfaceProperties& interface
)
{
    word alpharScheme("div(phirb,alpha)");
	word YiScheme("div(phi,Yi)");

    // Standard face-flux compression coefficient
    surfaceScalarField phic(mag(phi_ / U_.mesh().magSf()));

    surfaceScalarField phir(phic * interface.nHatf());

    // TODO: move this to general solve() function
    alpha2_ = 1 - alpha1_;

    // Soluables
    PtrList<volScalarField>& Si = reaction_->Y();
    // forAll(Si, i)
	// {
        label i = 6;

        volScalarField& Yi = Si[i];

        scalar maxYi = max(gMax(Yi), gMax(Yi.boundaryField())) + 1e-30;

        // normalizing
        Yi.oldTime() == Yi.oldTime() / maxYi;
        Yi == Yi / maxYi;

		surfaceScalarField phiComp = fvc::flux
        (
            -fvc::flux(-phir, alpha2_, alpharScheme),
            alpha1_,
            alpharScheme
        );

        tmp<surfaceScalarField> tYiPhi1Un
        (
            fvc::flux
            (
                phi_,
                Yi,
                YiScheme
            )
		  + phiComp * compressionCoeff(Yi)
        );

        {
            surfaceScalarField YiPhi10 = tYiPhi1Un;

            MULES::explicitSolve
            (
                geometricOneField(),
                Yi,
                phi_,
                YiPhi10,
                zeroField(),
                zeroField(),
                oneField(),
                zeroField()
            );
        }

        Yi.oldTime() == Yi.oldTime() * maxYi;
        Yi == Yi * maxYi;
    // }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admMixture::admMixture
(
    const volScalarField& Top,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    incompressibleTwoPhaseMixture(U, phi),
    H_
    (
        "test",
        dimless,
        this->subDict(get<wordList>("phases")[0]).lookupOrDefault
        (
            "H",
            1e-12
        )
    ),
    DS_
    (
        "test1",
        dimArea/dimTime,
        this->subDict(get<wordList>("phases")[0]).lookupOrDefault
        (
            "D",
            1e-8
        )
    ),
    DG_
    (
        "test2",
        dimArea/dimTime,
        this->subDict(get<wordList>("phases")[1]).lookupOrDefault
        (
            "D",
            1e-8
        )
    ),
    DalphaS_
    (
        IOobject
		(
			"DalphaS",
            alpha1_.time().timeName(),
			U_.mesh(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		U_.mesh(),
		dimensionedScalar
        (
            "DalphaSdefault",
            dimArea/dimTime,
            SMALL
        )
    ),
    DalphaG_
    (
        IOobject
		(
			"DalphaG",
            alpha1_.time().timeName(),
			U_.mesh(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		U_.mesh(),
		dimensionedScalar
        (
            "DalphaGdefault",
            dimArea/dimTime,
            SMALL
        )
    ),
    phiHS_
    (
        IOobject
		(
			"phiH",
            alpha1_.time().timeName(),
			U_.mesh(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		U_.mesh(),
		dimensionedScalar
        (
            "phiHdefault",
            dimVolume/dimTime,
            SMALL
        )
    ),
    reaction_
    (
        ADMno1::New(Top, U_.mesh())
    ),
    mDot_
    (
        IOobject
		(
			"mDot",
            alpha1_.time().timeName(),
			U_.mesh(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		U_.mesh(),
		dimensionedScalar
        (
            "mDotdefault",
            dimDensity/dimTime,
            SMALL
        )
    ),
    mDotAlphal_
    (
        IOobject
		(
			"mDotAlphal",
            alpha1_.time().timeName(),
			U_.mesh(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		U_.mesh(),
		dimensionedScalar
        (
            "mDotAlphaldefault",
            dimDensity/dimTime,
            SMALL
        )
    ),
    vDot_
    (
        IOobject
		(
			"vDot",
            alpha1_.time().timeName(),
			U_.mesh(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		U_.mesh(),
		dimensionedScalar
        (
            "vDotdefault",
            dimless/dimTime,
            SMALL
        )
    ),
    vDotAlphal_
    (
        IOobject
		(
			"vDotAlphal",
            alpha1_.time().timeName(),
			U_.mesh(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		U_.mesh(),
		dimensionedScalar
        (
            "vDotAlphaldefault",
            dimless/dimTime,
            SMALL
        )
    )
    // testing
    ,mDotTest_
    (
        "mDotTest",
        dimDensity/dimTime,
        -this->lookupOrDefault
        (
            "mDotTest",
            5e-3
        )
    ),
    phiD_
    (
        IOobject
        (
            "phiD",
            alpha1_.time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar
        (
            "phiD",
            dimMass/dimTime,
            SMALL
        )
    ),
    Mflux_
    (
        IOobject
        (
            "Mflux",
            alpha1_.time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
            // IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar
        (
            "MFlux",
            dimMass/dimVolume/dimTime,
            SMALL
        )
    )
{
    // alpha2_ = 1 - alpha1_;
    // alpha2_ = 1 - alpha1_;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::admMixture::limit()
{
    // Calculate mass transfer flux ----------------------------------------------
    // Info<< "correcting Yi" << endl;
    
    label i = 6;
    volScalarField& Yi = reaction_->Y()[i];

    //calculate alpha downwind
    surfaceScalarField fluxDir = fvc::snGrad(alpha1_)* U_.mesh().magSf();
    surfaceScalarField alphaDown = downwind<scalar>(U_.mesh(),fluxDir).interpolate(alpha1_);

    // Re-initialize transfer flux
    phiD_ = 0 * phiD_;

    // if (phiHScheme_ == "Gauss upwind")
    // {
    // 	forAll(species_, i)
    // 	{
    // 		phiD_+=Mw_[i]*
    // 			 (
    // 			 	DmY(i)*fvc::snGrad(Yi)*mesh.magSf()
    // 	           -fvc::flux(phiHUp(i),Yi,"div(phiHS,Yi)")
    // 	           -fvc::flux(phiHDown(i),Yi,"div(phiHS,Yi)")
    // 			 );
    // 	}
    // }
    // else if (phiHScheme_ == "Gauss linear")
    // {
    	// forAll(species_, i)
    	// {
            phiD_ += // Mw_*
            (
                DalphaS_ * fvc::snGrad(Yi) * U_.mesh().magSf()
              - fvc::flux(phiHS_, Yi, "div(phiH,Yi)")
            );
    	// }
    // }
    // else
    // {
    //     Info<< "div(phiHS,Yi) should be equal to Gauss linear or Gauss upwind"
    // 	<< endl
    // 	<< abort(FatalError);
    // }

    // Compute flux
    Mflux_ = fvc::div(phiD_ * alphaDown) - alpha1_ * fvc::div(phiD_);

    // alpha2_ = 1 - alpha1_;

    // compute Yi1 and Yi2
    // forAll(species_, i)
    // {
    //     volScalarField& Yi = Yi;
        // volScalarField& Y1i = phase1SpeciesMixture_.Y(i);
        // volScalarField& Y2i = phase2SpeciesMixture_.Y(i);
        // dimensionedScalar HYi = HY_[i];
        volScalarField Y1i = Yi / (alpha1_ + H_*(1 - alpha1_));
        volScalarField Y2i = H_ * Yi / (alpha1_ + H_*(1 - alpha1_));
    // }

    // // set saturation
    // phase1SpeciesMixture_.setSaturation(alpha1);
    // phase2SpeciesMixture_.setSaturation(alpha2);

    // // correct each phase
    // phase1SpeciesMixture_.correct();
    // phase2SpeciesMixture_.correct();

    // compute Yi from Y1i and Y2i
    // forAll(species_, i)
    // {
        // volScalarField& Yi = Yi;
        // volScalarField& Y1i = phase1SpeciesMixture_.Y(i);
        // volScalarField& Y2i = phase2SpeciesMixture_.Y(i);
        Yi = Y1i*alpha1_ + Y2i*(1 - alpha1_);
    // }

    // DEBUG
    Info<< "Yi correction: max(Mflux) = "
        << gMax(Mflux_.internalField())
        << "  Min(" << Yi.name() << ") = " << gMin(Yi.internalField())
        << "  Max(" << Yi.name() << ") = " << gMax(Yi.internalField())
        << "  Species concentration (sum) = "<< gSum(Yi.internalField())
        << endl;
};


const Foam::volScalarField&
Foam::admMixture::mDot()
{
    volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    mDot_ = limitedAlpha1*mDotTest_; 

    return mDot_;
}


const Foam::volScalarField&
Foam::admMixture::mDotAlphal()
{
    mDotAlphal_ = mDotTest_;
     
    return mDotAlphal_;
}


const Foam::volScalarField&
Foam::admMixture::vDot()
{
    dimensionedScalar pCoeff(1.0/this->rho1() - 1.0/this->rho2());

    vDot_ = pCoeff*this->mDot();

    return vDot_;
}


const Foam::volScalarField&
Foam::admMixture::vDotAlphal()
{
    volScalarField alphalCoeff
    (
        1.0/this->rho1() - this->alpha1()
       *(1.0/this->rho1() - 1.0/this->rho2())
    );

    vDotAlphal_ = alphalCoeff*this->mDotAlphal();

    return vDotAlphal_;
}


void Foam::admMixture::solve
(
    const interfaceProperties& interface
)
{
    // // DEBUG
    // Info<< ">>> testing admMixture::solve()" << endl;

    massTransferCoeffs();

    speciesMules(interface);
}


void Foam::admMixture::reaction
(
    const volScalarField& Top,
    const volScalarField& p
)
{
    reaction_->clear();

    reaction_->correct
    (
        this->phi_,
        this->alpha1_,
        this->alpha2_,
        Top,
        p
    );
}


// Foam::Pair<Foam::tmp<Foam::volScalarField>>
// Foam::admMixture::vDotAlphal() const
// {
//     volScalarField alphalCoeff
//     (
//         1.0/mixture_.rho1() - mixture_.alpha1()
//        *(1.0/mixture_.rho1() - 1.0/mixture_.rho2())
//     );

//     Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal();

//     return Pair<tmp<volScalarField>>
//     (
//         alphalCoeff*mDotAlphal[0],
//         alphalCoeff*mDotAlphal[1]
//     );
// }


// Foam::Pair<Foam::tmp<Foam::volScalarField>>
// Foam::admMixture::vDot() const
// {
//     dimensionedScalar pCoeff(1.0/mixture_.rho1() - 1.0/mixture_.rho2());
//     Pair<tmp<volScalarField>> mDot = this->mDot();

//     return Pair<tmp<volScalarField>>(pCoeff*mDot[0], pCoeff*mDot[1]);
// }


// bool Foam::admMixture::read()
// {
//     if (regIOobject::read())
//     {
//         return true;
//     }

//     return false;
// }


// ************************************************************************* //
