CorrectPhi
(
    U,
    phi,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU())),
    geometricZeroField(),
    pimple
);

#include "continuityErrs.H"

// CorrectPhi
// (
//     U,
//     phi,
//     p,
//     dimensionedScalar("rAUf", dimTime, 1),
//     geometricZeroField(),
//     pimple
// );


