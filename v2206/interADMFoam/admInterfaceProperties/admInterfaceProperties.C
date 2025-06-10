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

#include "admInterfaceProperties.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcAverage.H"
#include "fvCFD.H"
#include "unitConversion.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admInterfaceProperties::
admInterfaceProperties
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    interfaceProperties
    (
        alpha1,
        U,
        dict
    ),
    alpha1_(alpha1),
    //- Sharp force coefficient (put 0.98-0.99 for static problems, 0.4-0.5 for dynamic)
    cPc_
    (
        readScalar
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("cPc")
        )
    ),
    nNonOrthogonalCorrectors_
    (
        readLabel
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("nNonOrthogonalCorrectors")
        )
    ),
    pc_
    (
        IOobject
        (
            "pc",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh()
    ),
    //Specific cell where the capillary pressure value is used as refference
    pcRefCell_
    (
        readLabel
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("pcRefCell")
        )
    ),
    //value of the capillary pressure in the specific cell
    pcRefValue_
    (
        readScalar
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("pcRefValue")
        )
    ),
    phic_
    (
        IOobject
        (
            "phic",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("phic", dimPressure / dimLength*dimArea, 0.0)
    )
{

    setRefCell
    (
        pc_,
        alpha1_.mesh().solutionDict().subDict("PIMPLE"),
        pcRefCell_,
        pcRefValue_
    );

    interfaceProperties::correct();
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::admInterfaceProperties::calculatePhic()
{
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceScalarField& magSf = mesh.magSf();
 
    volScalarField alpha_pc
    (
        IOobject 
        (
            "alpha_pc",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_,
        pc_.boundaryField().types()
    );

    // Sharpen interface function
    // Raeini's thesis, equation (2.22) and (2.23)
    alpha_pc =
    (
        1.0 / (1.0 - cPc_) * (min(max(alpha1_, cPc_ / 2.0), (1.0 - cPc_ / 2.0)) - cPc_ / 2.0)
    );

    alpha_pc.correctBoundaryConditions();

    surfaceScalarField deltasf = fvc::snGrad(alpha_pc);

    //surface tension force
    surfaceScalarField stf = fvc::interpolate(sigmaK()) * deltasf;


    //surface tension force flux
    phic_ = stf*magSf;

    for(int nonOrth = 0; nonOrth <= nNonOrthogonalCorrectors_; nonOrth++)
    {
		//solve for pc
        fvScalarMatrix pcEqn
        (
            fvm::laplacian(pc_) == fvc::div(phic_)
        );

        pcEqn.setReference(pcRefCell_, pcRefValue_);

        solverPerformance residual = pcEqn.solve();

        
	    if (nonOrth==0) eqnResidual_ = residual.initialResidual();

		//add flux of pc to capillary flux
        if (nonOrth == nNonOrthogonalCorrectors_)
        { 
            phic_ -= pcEqn.flux();
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::admInterfaceProperties::correct()
{
    interfaceProperties::correct();

    calculatePhic();
}


// ************************************************************************* //

