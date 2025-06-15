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


// Create mixture and reaction
Info<< "Creating admMixture\n" << endl;
autoPtr<admMixture> mixture
(
    admMixture::New
    (
        Top, // <-- for ADMno1 reaction
        U,
        phi
    )
);

volScalarField& alpha1(mixture->alpha1());
volScalarField& alpha2(mixture->alpha2());

const dimensionedScalar& rho1 = mixture->rho1();
const dimensionedScalar& rho2 = mixture->rho2();


// Construct interface from alpha1 distribution
interfaceProperties interface
(
    alpha1,
    U,
    mixture()
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


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New
    (
        U, 
        phi,
        mixture()
    )
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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


// ---------------------------------------------------------------------------------
// TODO: implement when T is included
// Turbulent Prandtl number
// dimensionedScalar Prt("Prt", dimless, thermo->transportPropertiesDict());

// volScalarField kappaEff
// (
//     IOobject
//     (
//         "kappaEff",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     thermo->kappa()
// );

// Need to store rho for ddt(rhoCp, U)
// volScalarField rhoCp
// (
//     IOobject
//     (
//         "rhoCp",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     rho*thermo->Cp()
// );

// rhoCp.oldTime();

// volScalarField& T = thermo->T();
// ---------------------------------------------------------------------------------


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
PtrList<volScalarField>& Si = mixture->Si(); // <-- TODO: const?
PtrList<volScalarField>& Gi = mixture->Gi(); // <-- TODO: const?


