tmp<fvVectorMatrix> tUcEqn
(
    fvm::ddt(alphac, Uc) + fvm::div(alphaPhic, Uc)
  // - fvm::Sp(fvc::ddt(alphac) + fvc::div(alphaPhic), Uc)
  + continuousPhaseTurbulence->divDevRhoReff(Uc)
 ==
    (1.0/rhoc)*cloudSU
);
fvVectorMatrix& UcEqn = tUcEqn.ref();
UcEqn.relax();

volScalarField rAUc(1.0/UcEqn.A());
surfaceScalarField rAUcf("Dp", fvc::interpolate(rAUc));

surfaceScalarField phicForces
(
   fvc::flux(rAUc*cloudVolSUSu/rhoc) + rAUcf*(g & mesh.Sf())
);

if (pimple.momentumPredictor())
{
    solve
    (
        UcEqn
     ==
        fvc::reconstruct(phicForces/rAUcf)
      - fvc::grad(p)
    );
}


