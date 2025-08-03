{
	word alpharScheme("div(phirb,alpha)");
	word YiScheme("div(phi,Yi)");

    volScalarField& Yi = YPtrs[6];
    scalar HY = 0.1;

    // -------------------------------------------------------------------------------

    // {
        //Direction of interfacial flux
        surfaceScalarField fluxDir = fvc::snGrad(alphaLiq)*mesh.magSf();

        //upwind and downwind alpha1
        surfaceScalarField alphaUp = upwind<scalar>(mesh,fluxDir).interpolate(alphaLiq);
        surfaceScalarField alphaDown = downwind<scalar>(mesh,fluxDir).interpolate(alphaLiq);

        //upwind and downwnd Yi
        surfaceScalarField YiUp   = upwind<scalar>(mesh,fluxDir).interpolate(Yi);
        surfaceScalarField YiDown = downwind<scalar>(mesh,fluxDir).interpolate(Yi);
        
        dimensionedScalar sgn = sign(max(alphaDown*YiDown)-max((1-alphaUp)*YiUp));

        //normal compression coefficient
        surfaceScalarField deltaYi1 = max(-max(Yi),min(max(Yi),(YiDown-YiUp)/(alphaDown-alphaUp+1e-4)));
        
        //standard compression coefficient
        surfaceScalarField deltaYi2 = max(-max(Yi),min(max(Yi),YiDown/(alphaDown+(1-alphaDown)*HY)-HY*YiUp/(alphaUp+(1-alphaUp)*HY)));

        // return sgn*max(mag(deltaYi1),mag(deltaYi2));
        surfaceScalarField compressionCoeff = sgn*max(mag(deltaYi1),mag(deltaYi2));

    // }


    // -------------------------------------------------------------------------------

    // Standard face-flux compression coefficient
    surfaceScalarField phic(mag(phi/mesh.magSf()));

    // TODO: mixture.nHatf(alpha1_, alpha2_) instead of mixture.nHatf_) 
    surfaceScalarField phir(phic*mixture.nHatf_);

    // basicTwoPhaseMultiComponentTransportMixture& speciesMixture =pSpeciesMixture();
    // const speciesTable& solutionSpecies = speciesMixture.species();

    // forAll(solutionSpecies, i)
	// {
        // volScalarField& Yi = speciesMixture.Y(i);

        scalar maxYi = max(gMax(Yi),gMax(Yi.boundaryField()))+1e-30;

        Yi.oldTime() == Yi.oldTime()/maxYi;
        Yi == Yi/maxYi;

		surfaceScalarField phiComp = fvc::flux
        (
            -fvc::flux(-phir, alphaGas, alpharScheme),
            alphaLiq,
            alpharScheme
        );

        tmp<surfaceScalarField> tYiPhi1Un
        (
            fvc::flux
            (
                phi,
                Yi,
                YiScheme
            )
		 +  phiComp*compressionCoeff //*speciesMixture.compressionCoeff(i)
        );

        {
            surfaceScalarField YiPhi10 = tYiPhi1Un;

            MULES::explicitSolve
            (
                geometricOneField(),
                Yi,
                phi,
                YiPhi10,
                zeroField(),
                zeroField(),
                oneField(),
                zeroField()
            );
        }

        Yi.oldTime() == Yi.oldTime()*maxYi;
        Yi == Yi*maxYi;
    // }
}
