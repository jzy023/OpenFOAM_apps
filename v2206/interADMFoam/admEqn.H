Info<< "ADM convection-diffusion calculation" << endl;

// TODO: implement this with fields tabble in createFields.H
//       to enable divSchemes for Si and Gi [ref: reactingFoam]

// tmp<fv::convectionScheme<scalar>> mvConvection
// (
//     fv::convectionScheme<scalar>::New
//     (
//         mesh,
//         fields,
//         phi,
//         mesh.divScheme("div(phi,Yi_h)")
//     )
// );


//TODO: turbulent diff and different diff compared to particulate species
volScalarField Deff = 
(
    (turbulence->nu() / reaction->Sc())
  + (turbulence->nut() / reaction->Sct())
);

{
    volScalarField& Yi = reaction->Y()[6];
    // volScalarField& Yi = YPtrs[6];
    // volScalarField& Yi = speciesMixture.Y(i);

    //compute diffusion coefficient
    // surfaceScalarField D = speciesMixture.DmY(i);

    //calculate Henry's transfer flux
    // surfaceScalarField phiH = speciesMixture.phiH(i);
    // surfaceScalarField phiHUp = speciesMixture.phiHUp(i);
    // surfaceScalarField phiHDown = speciesMixture.phiHDown(i);
    
    // surfaceScalarField phiH = DmY(i)*(1-HY_[i])/fvc::interpolate((alpha1_+HY_[i]*(1-alpha1_)))*fvc::snGrad(alpha1_)*mesh_.magSf();
    dimensionedScalar DmY
    (
        Deff.dimensions(),
        1e-16
    );
    scalar HY = 0.1;
    surfaceScalarField phiH = 
    (
        DmY * (1 - HY) / fvc::interpolate((alpha1 + HY * (1 - alpha1)))
      * fvc::snGrad(alpha1)*mesh.magSf()
    );

    volScalarField SpY = (Yi-Yi.oldTime())/runTime.deltaT();

    while (pimple.correctNonOrthogonal())
    {
        //compute Y equation
        fvScalarMatrix YiEqn
        (
            fvm::ddt(Yi)
          + fvm::div(phiH,Yi,"div(phiH,Yi)")
          - fvm::laplacian(DmY,Yi)
        //   - fvm::laplacian(Deff,Yi)
          ==
            SpY
        );

        // if (phiHScheme == "Gauss linear") YiEqn += fvm::div(phiH,Yi,"div(phiH,Yi)");
        // else if (phiHScheme == "Gauss upwind")
        // {
        //     YiEqn += fvm::div(phiHUp,Yi,"div(phiH,Yi)")
        //     + fvm::div(phiHDown,Yi,"div(phiH,Yi)");
        // }
        // else
        // {
        //     Info<< "div(phiH,Yi) should be equal to Gauss linear or Gauss upwind"
		// 	<< endl
		// 	<< abort(FatalError);
        // }

        //solve equations
        if (pimple.finalNonOrthogonalIter())
        {
            YiEqn.solve(mesh.solver("YiFinal"));
        }
        else
        {
            YiEqn.solve(mesh.solver("Yi"));
        }
    }

    // Info<< "Species concentration = "
    // << Yi.weightedAverage(mesh.V()).value()
    Info<< "Species concentration (sum) = "
    << gSum(Yi.internalField())
    << "  Min(" << Yi.name() << ") = " << gMin(Yi.internalField())
    << "  Max(" << Yi.name() << ") = " << gMax(Yi.internalField())
    << endl;

    //calculate mass transfer flux ----------------------------------------------
    Info<< "correcting Yi" << endl;
    // speciesMixture.correct();

    //calculate alpha downwind
    surfaceScalarField fluxDir = fvc::snGrad(alpha1)*mesh.magSf();
    surfaceScalarField alphaDown = downwind<scalar>(mesh,fluxDir).interpolate(alpha1);

    //Re-initialize transfer flux
	phiD_ =0*phiD_;

    // if (phiHScheme_ == "Gauss upwind")
    // {
	// 	forAll(species_, i)
	// 	{
	// 		phiD_+=Mw_[i]*
	// 			 (
	// 			 	DmY(i)*fvc::snGrad(Yi)*mesh.magSf()
	// 	           -fvc::flux(phiHUp(i),Yi,"div(phiH,Yi)")
	// 	           -fvc::flux(phiHDown(i),Yi,"div(phiH,Yi)")
	// 			 );
	// 	}
    // }
    // else if (phiHScheme_ == "Gauss linear")
    // {
	// 	forAll(species_, i)
	// 	{
			phiD_+= // Mw_*
				 (
				 	DmY*fvc::snGrad(Yi)*mesh.magSf()
				   -fvc::flux(phiH,Yi,"div(phiH,Yi)")
				 );
	// 	}
    // }
    // else
    // {
    //     Info<< "div(phiH,Yi) should be equal to Gauss linear or Gauss upwind"
	// 	<< endl
	// 	<< abort(FatalError);
    // }

    //Compute flux
	Mflux_ = fvc::div(phiD_*alphaDown)-alpha1*fvc::div(phiD_);

    volScalarField alpha2 = 1-alpha1;

	//compute Yi1 and Yi2
	// forAll(species_, i)
	// {
		// volScalarField& Yi = Yi;
		// volScalarField& Y1i = phase1SpeciesMixture_.Y(i);
		// volScalarField& Y2i = phase2SpeciesMixture_.Y(i);
		// dimensionedScalar HYi = HY_[i];
		Y1i = Yi/(alpha1+HY*(1-alpha1));
		Y2i = HY*Yi/(alpha1+HY*(1-alpha1));
	// }

    // //set saturation
    // phase1SpeciesMixture_.setSaturation(alpha1);
    // phase2SpeciesMixture_.setSaturation(alpha2);

    // //correct each phase
    // phase1SpeciesMixture_.correct();
    // phase2SpeciesMixture_.correct();

	//compute Yi from Y1i and Y2i
	// forAll(species_, i)
	// {
		// volScalarField& Yi = Yi;
		// volScalarField& Y1i = phase1SpeciesMixture_.Y(i);
		// volScalarField& Y2i = phase2SpeciesMixture_.Y(i);
		Yi = Y1i*alpha1 + Y2i*(1-alpha1);
	// }

    // DEBUG
    Info<< "max(Mflux) = "
        << gMax(Mflux_.internalField())
        << endl;

    // --------------------------------------------------------------------------------------------------

    // volScalarField& Yi = YPtrs[6];
    // for (int i = 0; i < 3; i++)
    // {
    //     fvScalarMatrix YiEqn
    //     (
    //         fvm::ddt(Yi)
    //       + fvm::div(phi, Yi, "div(phi,test)")
    //       - fvm::laplacian(Deff, Yi)
    //     ==
    //         alpha1*fvOptions(Yi)
    //     );

    //     YiEqn.relax();
    //     fvOptions.constrain(YiEqn);
    //     YiEqn.solve();
    //     fvOptions.correct(Yi);

    //     tPhiYi = YiEqn.flux();
    // }

    // MULES::explicitSolve
    // (
    //     geometricOneField(),
    //     Yi,
    //     phi,
    //     tPhiYi.ref(),
    //     oneField(),
    //     zeroField()
    // );

    // --------------------------------------------------------------------------------------------------

        // volScalarField& Yi = YPtrs[6];

        // // Deff *= pos(alpha1 - 0.99);
        // // tmp<surfaceScalarField> tPhiYi;

        // for (int i = 0; i < 2; i++)
        // {
        //     fvScalarMatrix YiEqn
        //     (
        //         fvm::ddt(Yi)
        //       - fvm::laplacian(alpha1*Deff, Yi)
        //       + fvm::div(mixture.phases()["liquid"].alphaPhi(), Yi, "div(phi,test)")
        //     //   + fvm::div(phi, Yi, "div(phi,test)")
        //     ==
        //         alpha1*fvOptions(Yi)
        //     );

        //     YiEqn.relax();
        //     fvOptions.constrain(YiEqn);
        //     YiEqn.solve();
        //     fvOptions.correct(Yi);

        //     // tPhiYi = YiEqn.flux();
        // }

        // MULES::explicitSolve
        // (
        //     geometricOneField(),
        //     Yi,
        //     phi,
        //     tPhiYi.ref(),
        //     oneField(),
        //     zeroField()
        // );

    // --------------------------------------------------------------------------------------------------

        //- Fluid and Particulate Phase

        // const volScalarField residualAlpha_
        // (
        //     IOobject
        //     (
        //         "residualAlpha",
        //         mesh.time().timeName(),
        //         mesh,
        //         IOobject::NO_READ, // READ_IF_PRESENT,
        //         IOobject::NO_WRITE
        //     ),
        //     mesh,
        //     dimensionedScalar
        //     (
        //         dimless,
        //         1e-3
        //     )
        // );
            
        // volScalarField& Yi = YPtrs[6];

        // // Deff *= pos(alpha1 - 0.99);

        // for (size_t i = 0; i < 3; i++)
        // {
        //     fvScalarMatrix YiEqn
        //     (
        //         fvm::ddt(alpha1, Yi)
        //       - fvm::laplacian(alpha1*Deff, Yi)
        //       + fvm::div(mixture.phases()["liquid"].alphaPhi(), Yi, "div(phi,test)")
        //     //   + fvm::div(phi, Yi, "div(phi,test)")
        //     ==
        //     //     alpha1*fvOptions(Yi)
        //        fvc::ddt(residualAlpha_, Yi)
        //      - fvm::ddt(residualAlpha_, Yi)
        //     );

        //     YiEqn.relax();
        //     fvOptions.constrain(YiEqn);
        //     YiEqn.solve();
        //     fvOptions.correct(Yi);
        // }


    // --------------------------------------------------------------------------------------------------

        // - Fluid and Particulate Phase
        
        // volScalarField& Yi = YPtrs[6];
        // volScalarField YiAlpha = Yi;
        // YiAlpha.field() = max(alpha1, scalar(1e-3))().field() * Yi.field();

        // Deff *= pos(alpha1 - 0.99);

        // for (size_t i = 0; i < 3; i++)
        // {
        //     fvScalarMatrix YiEqn
        //     (
        //         fvm::ddt(YiAlpha)
        //     //   + fvm::div(phi, YiAlpha, "div(phi,test)")
        //       + fvm::div(mixture.phases()["liquid"].alphaPhi(), Yi, "div(phi,test)")
        //       - fvm::laplacian(Deff, YiAlpha)
        //     ==
        //         alpha1*fvOptions(YiAlpha)
        //     );

        //     YiEqn.relax();
        //     fvOptions.constrain(YiEqn);
        //     YiEqn.solve();
        //     fvOptions.correct(YiAlpha);
        //     Yi = YiAlpha / max(alpha1, scalar(1e-3));
        // }
    

    // --------------------------------------------------------------------------------------------------

    // const volScalarField residualAlpha_
    // (
    //     IOobject
    //     (
    //         "residualAlpha",
    //         mesh.time().timeName(),
    //         mesh,
    //         IOobject::NO_READ, // READ_IF_PRESENT,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar
    //     (
    //         dimless,
    //         1e-3
    //     )
    // );

    // const volScalarField rhoLiq
    // (
    //     IOobject
    //     (
    //         "rhoLiq",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     mixture.phases()["liquid"].rho()
    // );

    // const surfaceScalarField& alphaRhoPhiLiq = fvc::interpolate(rhoLiq)*mixture.phases()["liquid"].alphaPhi();

    // //- Fluid and Particulate Phase
        
    //     volScalarField& Yi = YPtrs[6];

    //     fvScalarMatrix YiEqn
    //     (
    //         fvm::ddt(alpha1, rhoLiq, Yi)
    //       + fvm::div(alphaRhoPhiLiq, Yi, "div(alphaRhoPhi.*,Yi)")
    //       - fvm::laplacian
    //         (
    //             fvc::interpolate(alpha1)*
    //             fvc::interpolate(rhoLiq*Deff),
    //             Yi
    //         )
    //      ==
    //         fvc::ddt(residualAlpha_*rhoLiq, Yi)
    //       - fvm::ddt(residualAlpha_*rhoLiq, Yi)
    //     );

    //     YiEqn.relax();
    //     fvOptions.constrain(YiEqn);
    //     YiEqn.solve();
    //     fvOptions.correct(Yi);

}