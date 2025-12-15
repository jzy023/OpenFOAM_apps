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

//- Find isCellsInterface_
void Foam::admMixture::kLaCells()
{
    isCellsInterface_ = 
    (
        max(zeroField(), alpha1_ - (1.0 - alphaInterface_))/(alpha1_ - (1.0 - alphaInterface_))
      + max(zeroField(), alphaInterface_ - alpha1_)/(alphaInterface_ - alpha1_)
    );

    isCellsInterface_ = 1.0 - isCellsInterface_;

    volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );
    
    // TODO: fix kLa for benchmark case?
    kLaCells_.field() = // limitedAlpha1 * // <- check this?
    (
        alphaI_*isCellsInterface_.field() + alphaW_*isCellsActWall_.field()
        // alphaI*isCellsInterface_.field() + (1/alphaW_.value())*isCellsActWall_.field()
    );
}

        
//- Find isCellsActWall_
void Foam::admMixture::findCellsActWall()
{
    forAll(actPatch_, wallI)
    {
        // To access the boundary patches information
        const fvPatch& cPatch = U_.mesh().boundary()[actPatch_[wallI]];

        // Starting index of the face in a patch
        label faceId_start = cPatch.start() ;

        // List of cells close to a boundary
        const labelUList& faceCells = cPatch.faceCells();

        forAll(cPatch, faceI) 
        { 
            // index of each face
            label faceID = faceId_start + faceI;

            // id of the owner cell having the face
            label faceOwner = faceCells[faceI];

            // mark the cells
            isCellsActWall_[faceOwner] = 1;
        }
    }
}


//- Species transport equations
void Foam::admMixture::massTransferCoeffs()
{
    volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    volScalarField limitedAlpha2
    (
        min(max(alpha2_, scalar(0)), scalar(1))
    );
            
    // calculate and return mean diffusion coefficient
    // TODO: multispecies? add turbulent diffusivity too?
    D2Eff_ = fvc::interpolate(D2_ * limitedAlpha2);
    D1Eff_ = 
    (
        fvc::interpolate(D1_ * limitedAlpha1)
        // fvc::interpolate(D1_ * limitedAlpha1 + H_ * DS2_ * limitedAlpha2) / fvc::interpolate(limitedAlpha1 + H_ * limitedAlpha2)
    );
            
    //- Calculate interface mass transfer flux by Henry's Law
    // TODO: check for gas [all 0 if H_ -> 0]
    // surfaceScalarField phiHUp = speciesMixture.phiHUp(i);
    // surfaceScalarField phiHDown = speciesMixture.phiHDown(i);
    phiHS_ = 
    (
        D1Eff_ * (1 - H_) / fvc::interpolate((limitedAlpha1 + H_ * (1 - limitedAlpha1)))
      * fvc::snGrad(limitedAlpha1) * U_.mesh().magSf()
    );
}


//- interface compreassion coefficient
surfaceScalarField Foam::admMixture::compressionCoeff
(
    const volScalarField& Yi
)
{
    volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    // Direction of interfacial flux
    surfaceScalarField fluxDir = fvc::snGrad(limitedAlpha1) * U_.mesh().magSf();

    // Upwind and downwind alpha1
    surfaceScalarField alphaUp = upwind<scalar>(U_.mesh(),fluxDir).interpolate(limitedAlpha1);
    surfaceScalarField alphaDown = downwind<scalar>(U_.mesh(),fluxDir).interpolate(limitedAlpha1);

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
    
    // Soluables
    forAll(SiAlpha_, i)
	{
        if (i == 7)
        {
            // divPhi = fvc::div(phiHS_);
            // divPhiSh2 = fvc::div(phiHS_, Sh2);
            continue;
        }

        volScalarField& Yi = SiAlpha_[i];
        
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
    }
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
    reaction_
    (
        ADMno1::New
        (
            this->alpha1_,
            Top, 
            U_.mesh()
        )
    ),
    // alphaW_
    // (
    //     "alphaW",
    //     dimless,
    //     this->subDict("degassing").get<scalar>("alphaW")
    // ),
    alphaI_
    (
        this->subDict("degassing").lookupOrDefault
        (
            "alphaI",
            0.1
        )
    ),
    alphaW_
    (
        this->subDict("degassing").lookupOrDefault
        (
            "alphaW",
            0.12
        )
    ),
    R_
    (
        "R",
        dimPressure*dimVolume/dimMass/dimTemperature,
        this->subDict("degassing").lookupOrDefault
        (
            "R",
            8.3145
        )
    ),
    H_
    (
        "H",
        dimless,
        this->subDict("degassing").lookupOrDefault
        (
            "H",
            1e-12
        )
    ),
    D1_
    (
        "D1",
        dimArea/dimTime,
        this->subDict(get<wordList>("phases")[0]).lookupOrDefault
        (
            "D",
            1e-8
        )
    ),
    D2_
    (
        "D2",
        dimArea/dimTime,
        this->subDict(get<wordList>("phases")[1]).lookupOrDefault
        (
            "D",
            1e-8
        )
    ),
    D1Eff_
    (
        IOobject
		(
			"DSEff",
            alpha1_.time().timeName(),
			U_.mesh(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		U_.mesh(),
		dimensionedScalar
        (
            dimArea/dimTime,
            SMALL
        )
    ),
    D2Eff_
    (
        IOobject
		(
			"DGEff",
            alpha1_.time().timeName(),
			U_.mesh(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		U_.mesh(),
		dimensionedScalar
        (
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
            dimVolume/dimTime,
            SMALL
        )
    ),
    alphaInterface_
    (
        "interfaceThreshold",
        dimless,
        this->subDict("degassing").lookupOrDefault
        (
            "alphaInterface",
            0.8
        )
    ),
    actPatch_
    (
        this->subDict("degassing").subDict("walls").toc()
    ),
    isCellsFull_
    (
        IOobject
        (
            "isCellsFull_",
            alpha1_.time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar
        (
            dimless,
            Zero
        )
    ),
    // testing -------------------------------------------------------
    isCellsEmpty_
    (
        IOobject
        (
            "isCellsEmpty_",
            alpha1_.time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar
        (
            dimless,
            Zero
        )
    ),
    isCellsInterface_
    (
        IOobject
        (
            "isCellsInterface_",
            alpha1_.time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar
        (
            dimless,
            Zero
        )
    ),
    // ---------------------------------------------------------------
    isCellsActWall_
    (
        IOobject
        (
            "isActWallCells",
            alpha1_.time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar
        (
            dimless,
            Zero
        )
    ),
    kLaCells_
    (
        IOobject
        (
            "kLaCells_",
            alpha1_.time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar
        (
            dimless/dimTime,
            Zero
        )
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
            Zero
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
            Zero
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
            Zero
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
            Zero
        )
    )
    // testing
    ,mDotTest_
    (
        "mDotTest",
        // dimDensity/dimTime,
        dimDensity,
        this->subDict("degassing").lookupOrDefault
        (
            "mDotTest",
            1e-12
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
    // ,divPhi
    // (
    //     IOobject
    //     (
    //         "divPhi",
    //         alpha1_.time().timeName(),
    //         U_.mesh(),
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     U_.mesh(),
    //     dimensionedScalar
    //     (
    //         "MFlux",
    //         dimless/dimTime,
    //         Zero
    //     )
    // ),
    // divPhiSh2
    // (
    //     IOobject
    //     (
    //         "divPhiSh2",
    //         alpha1_.time().timeName(),
    //         U_.mesh(),
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     U_.mesh(),
    //     dimensionedScalar
    //     (
    //         "MFlux",
    //         dimDensity/dimTime,
    //         Zero
    //     )
    // )
{
    //- Main substances concentration initialization

    // Info<< "Reading ADM no1 initial concentrations for soluables" << endl;
    
    label iNames = 0;

    SiAlpha_.resize(reaction_->Y().size());

    forAll(reaction_->namesSoluable, i)
    {
        SiAlpha_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    reaction_->namesSoluable[i],
                    alpha1_.time().timeName(),
                    U_.mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                U_.mesh()
            )
        );
    }

    iNames += reaction_->namesSoluable.size();

    forAll(reaction_->namesParticulate, i)
    {
        SiAlpha_.set
        (
            i + iNames,
            new volScalarField
            (
                IOobject
                (
                    reaction_->namesParticulate[i],
                    alpha1_.time().timeName(),
                    U_.mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                U_.mesh()
            )
        );
    }

    //- Gaseuoses initialization

    // Info<< "Reading ADM no1 initial concentrations for gaseuoses" << endl;

    GiAlpha_.resize(reaction_->GAve().size());

    forAll(reaction_->namesGaseous, i)
    {
        GiAlpha_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    reaction_->namesGaseous[i],
                    alpha1_.time().timeName(),
                    U_.mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                U_.mesh()
            )
        );
    }

    // testing
    limitAlpha();

    // initializing Si and Gi for ADMno1 reaction
    speciesADMCorrect();

    // initializing ADMno1 reaction
    reaction_->init(Top);

    // marking inter-phase mass transfer surfaces
    findCellsActWall();

    // initializing inter phase mass transfer rate [s-1]
    kLaCells();

    // // DEBUG
    // Info<< "H: "    << H_ 
    //     << "\nDS: " << D1_
    //     << "\nDG: " << D2_ 
    //     << "\nmDotTest: " << mDotTest_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Use this to log phase volume corrected generation rate
// TODO: make for Gh2, Gco2 and Gch4
const Foam::volScalarField&
Foam::admMixture::mDot()
{
    volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    // // ----------------------------------------------------------------------------------
    // mDot_ = limitedAlpha1 * 
    // (
    // //   - reaction_->GRAve()[0] - reaction_->GRAve()[1]
    // //   -(reaction_->GRAve()[2] * 44 / 1000) 
    //   - reaction_->GR()[0] - reaction_->GR()[1]
    //   -(reaction_->GR()[2] * 44 / 1000) 
    // );

    // // DEBUG
    // // volScalarField mDotControlled = limitedAlpha1 * mDotTest_ * kLaCells_;

    // // Info<< ">>> total gas generation rate [mol * m-3]: " << mDot_.weightedAverage(limitedAlpha1.mesh().V()).value() << " , "<< endl;
    // // Info<< ">>> test gas generation rate [mol * m-3]: " << mDotControlled.weightedAverage(limitedAlpha1.mesh().V()).value() << endl;     

    // return mDot_;

    // ----------------------------------------------------------------------------------
    // !!! Case studying with user-forced mDotTest_
    mDot_ = limitedAlpha1 * mDotTest_ * kLaCells_;

    return mDot_;
}


// use this for updateing dY in reaction_
// TODO:  make for Gh2, Gco2 and Gch4
const Foam::volScalarField&
Foam::admMixture::mDotAlphal()
{
    // // ----------------------------------------------------------------------------------
    // mDotAlphal_ =
    // (
    // //   - reaction_->GRAve()[0] - reaction_->GRAve()[1]
    // //   -(reaction_->GRAve()[2] * 44 / 1000) 
    //   - reaction_->GR()[0] - reaction_->GR()[1]
    //   -(reaction_->GR()[2] * 44 / 1000) 
    // );

    // return mDotAlphal_;

    // ----------------------------------------------------------------------------------
    // !!! Case studying with user-forced mDotAlphal_
    mDotAlphal_ = mDotTest_ * kLaCells_;

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


// testing
void Foam::admMixture::limitAlpha()
{
    forAll(this->alpha1(), celli)
    {
        if (this->alpha1()[celli] <= 0.0)
        {
            this->alpha1()[celli] = 1e-16;
        }
        else if (this->alpha1()[celli] >= 1.0)
        {
            this->alpha1()[celli] = 1 - 1e-16;
        }
    }
    this->alpha1().correctBoundaryConditions();

    this->alpha2() = 1 - this->alpha1();
    this->alpha2().correctBoundaryConditions();
}


void Foam::admMixture::solvePhase
(
    const interfaceProperties& interface
)
{
    // volScalarField limitedAlpha1
    // (
    //     min(max(alpha1_, scalar(0)), scalar(1))
    // );
    // SiAlpha_[7] = limitedAlpha1*reaction_->Y()[7];
    // SiAlpha_[7].correctBoundaryConditions();

    // scalar maxYiInternal = gMax(SiAlpha_[7]) + 1e-30;
    // scalar maxYiAll = max(gMax(SiAlpha_[7]), gMax(SiAlpha_[7].boundaryField())) + 1e-30;
    // Info<< ">>> Sh2 internal = " << maxYiInternal 
    //     << " , Sh2 all = " << maxYiAll << endl;

    alpha2_ = 1 - alpha1_;

    // Classify cells
    kLaCells();
    
    // Calculate interface mass flux: phiHS_
    massTransferCoeffs();

    // Solve MULES equations for each Yi
    speciesMules(interface);

    // Convert mesh based Yi to phase based Yi for ADM calculation
    speciesADMCorrect();

    // printGasGenRate();
}


void Foam::admMixture::solveReaction
(
    const surfaceScalarField& phi,
    const volScalarField& Top
)
{
    // // DEBUG
    // Info<< ">>> " << reaction_->PgasAve() << endl;
    
    reaction_->clear();

    reaction_->correct
    (
        kLaCells_,
        phi,
        Top
    );
}


void Foam::admMixture::updateReactionGas
(
    const dimensionedScalar deltaT
)
{
    PtrList<dimensionedScalar>& G = reaction_->GAve();
    PtrList<dimensionedScalar>& dG = reaction_->dGAve();
    forAll(dG, i)
    {
        G[i] += dG[i] * deltaT;
        GiAlpha_[i] = this->alpha2() * G[i];

        // DEBUG
        // Info<< G[i].name()  << " concentration = "
        //     << G[i].value() << ", generation rate =  "
        //     << reaction_->GRAve().value() << endl;
    }
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


// void Foam::admMixture::limit()
// {
//     // Calculate mass transfer flux ----------------------------------------------
//     // Info<< "correcting Yi" << endl;
//     // TODO:
//     label i = 6;
//     volScalarField& Yi = reaction_->Y()[i];

//     //calculate alpha downwind
//     surfaceScalarField fluxDir = fvc::snGrad(alpha1_)* U_.mesh().magSf();
//     surfaceScalarField alphaDown = downwind<scalar>(U_.mesh(),fluxDir).interpolate(alpha1_);

//     // Re-initialize transfer flux
//     phiD_ = 0 * phiD_;

//     // if (phiHScheme_ == "Gauss upwind")
//     // {
//     // 	forAll(species_, i)
//     // 	{
//     // 		phiD_+=Mw_[i]*
//     // 			 (
//     // 			 	DmY(i)*fvc::snGrad(Yi)*mesh.magSf()
//     // 	           -fvc::flux(phiHUp(i),Yi,"div(phiHS,Yi)")
//     // 	           -fvc::flux(phiHDown(i),Yi,"div(phiHS,Yi)")
//     // 			 );
//     // 	}
//     // }
//     // else if (phiHScheme_ == "Gauss linear")
//     // {
//     	// forAll(species_, i)
//     	// {
//             phiD_ += // Mw_*
//             (
//                 D1Eff_ * fvc::snGrad(Yi) * U_.mesh().magSf()
//               - fvc::flux(phiHS_, Yi, "div(phiH,Yi)")
//             );
//     	// }
//     // }
//     // else
//     // {
//     //     Info<< "div(phiHS,Yi) should be equal to Gauss linear or Gauss upwind"
//     // 	<< endl
//     // 	<< abort(FatalError);
//     // }

//     // Compute flux
//     Mflux_ = fvc::div(phiD_ * alphaDown) - alpha1_ * fvc::div(phiD_);

//     // alpha2_ = 1 - alpha1_;

//     // compute Yi1 and Yi2
//     // forAll(species_, i)
//     // {
//     //     volScalarField& Yi = Yi;
//         // volScalarField& Y1i = phase1SpeciesMixture_.Y(i);
//         // volScalarField& Y2i = phase2SpeciesMixture_.Y(i);
//         // dimensionedScalar HYi = HY_[i];
//         volScalarField Y1i = Yi / (alpha1_ + H_*(1 - alpha1_));
//         volScalarField Y2i = H_ * Yi / (alpha1_ + H_*(1 - alpha1_));
//     // }

//     // // set saturation
//     // phase1SpeciesMixture_.setSaturation(alpha1);
//     // phase2SpeciesMixture_.setSaturation(alpha2);

//     // // correct each phase
//     // phase1SpeciesMixture_.correct();
//     // phase2SpeciesMixture_.correct();

//     // compute Yi from Y1i and Y2i
//     // forAll(species_, i)
//     // {
//         // volScalarField& Yi = Yi;
//         // volScalarField& Y1i = phase1SpeciesMixture_.Y(i);
//         // volScalarField& Y2i = phase2SpeciesMixture_.Y(i);
//         Yi = Y1i*alpha1_ + Y2i*(1 - alpha1_);
//     // }

//     // // DEBUG
//     // Info<< "Yi correction: max(Mflux) = "
//     //     << gMax(Mflux_.internalField())
//     //     << "  Min(" << Yi.name() << ") = " << gMin(Yi.internalField())
//     //     << "  Max(" << Yi.name() << ") = " << gMax(Yi.internalField())
//     //     << "  Species concentration (sum) = "<< gSum(Yi.internalField())
//     //     << endl;
// };

// ************************************************************************* //
