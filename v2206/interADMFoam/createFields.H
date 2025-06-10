#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


// ADM reactions
Info<< "Reading field operational temperature" << endl;
volScalarField Top
(
    IOobject
    (
        "Top",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "TopDefault", 
        dimTemperature,
        308.15
    )
);

autoPtr<ADMno1> reaction
(
    ADMno1::New(Top, mesh)
);


// TODO: integrate and cleanup 
// Creating e based thermo
autoPtr<twoPhaseMixtureEThermo> thermo
(
    new twoPhaseMixtureEThermo(U, phi)
);

volScalarField& alpha1(thermo->alpha1());
volScalarField& alpha2(thermo->alpha2());

const dimensionedScalar& rho1 = thermo->rho1();
const dimensionedScalar& rho2 = thermo->rho2();


// Create mixture
Info<< "Creating admMixture\n" << endl;
autoPtr<admMixture> mixture
(
    admMixture::New
    (
        mesh,
        alpha1,
        reaction->Y(),
        reaction->G(),
        phi
    )
);


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

// Construct interface from alpha1 distribution
interfaceProperties interface
(
    alpha1,
    U,
    // TODO: fix this dict too
    thermo->transportPropertiesDict()
);


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, thermo())
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField& p = thermo->p();

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, thermo->transportPropertiesDict());

volScalarField kappaEff
(
    IOobject
    (
        "kappaEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo->kappa()
);

// Need to store rho for ddt(rhoCp, U)
volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*thermo->Cp()
);

rhoCp.oldTime();


// testing
surfaceScalarField phiD_
(
	IOobject
	(
		"phiD",
        runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("phiD", dimMass/dimTime, 0.0)
);

volScalarField Mflux_
(
	IOobject
	(
		"Mflux",
        runTime.timeName(),
		mesh,
		IOobject::NO_READ,
        IOobject::NO_WRITE
		// IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("MFlux", dimMass/dimVolume/dimTime, 0.0)
);

// 
PtrList<volScalarField>& Si = reaction->Y(); // <-- TODO: const?
PtrList<volScalarField>& Gi = reaction->G(); // <-- TODO: const?


