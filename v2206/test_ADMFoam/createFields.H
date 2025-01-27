#include "createRDeltaT.H"

// ADM reactions
Info<< "Reading field operational temperature" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "TDefault", 
        dimTemperature,
        308.15
    )
);

autoPtr<ADMno1> reaction
(
    ADMno1::New(T, mesh)
);


// Original [buoyantBoussinesqPimpleFoam]
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

// #include "readTransportProperties.H"
// singlePhaseTransportModel laminarTransport(U, phi);
// Info<< "Creating turbulence model\n" << endl;
// autoPtr<incompressible::turbulenceModel> turbulence
// (
//     incompressible::turbulenceModel::New(U, phi, laminarTransport)
// );

// // Kinematic density for buoyancy force
// volScalarField rhok
// (
//     IOobject
//     (
//         "rhok",
//         runTime.timeName(),
//         mesh
//     ),
//     1.0 - beta*(T - TRef)
// );

// Kinematic turbulent thermal thermal conductivity m2/s
Info<< "Reading field alphat\n" << endl;
volScalarField alphat
(
    IOobject
    (
        "alphat",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

// // Pressure field
// volScalarField p
// (
//     IOobject
//     (
//         "p",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     // p_rgh + rhok*gh
//     p_rgh + rho*gh
// );

// label pRefCell = 0;
// scalar pRefValue = 0.0;
// setRefCell
// (
//     p,
//     p_rgh,
//     pimple.dict(),
//     pRefCell,
//     pRefValue
// );

// if (p_rgh.needReference())
// {
//     p += dimensionedScalar
//     (
//         "p",
//         p.dimensions(),
//         pRefValue - getRefCellValue(p, pRefCell)
//     );
// }

// mesh.setFluxRequired(p_rgh.name());

#include "createMRF.H"
#include "createIncompressibleRadiationModel.H"
#include "createFvOptions.H"


// testing
// > multiphaseInterFoam
multiphaseMixture mixture(U, phi);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE 
    ),
    mixture.rho()
);
rho.oldTime();

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);

// Pressure field
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
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
}

mesh.setFluxRequired(p_rgh.name());






