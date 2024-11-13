/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |  
     \\/     M anipulation  |
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

#include "ADMno1.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ADMno1::gasPressure()
{
    //- gas pressure
    // volScalarField Ph2o = GPtrs_[0];

    volScalarField Ph2o
    (
        IOobject
        (
            "Ph2o",
            fac_.mesh().time().timeName(),
            fac_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fac_.mesh(),
        dimensionedScalar
        (
           "Ph2oDefault", 
            GPtrs_[0].dimensions(),
            Zero
        )
    );

    Ph2o.field() = para_.KH().h2o * exp(5290.0 * fac_ * 100 * R_);
    Pgas_.field() = (GPtrs_[0] / 16.0 + GPtrs_[1] / 64.0 + GPtrs_[2]) * R_ * TopDummy_ + Ph2o;
}


void Foam::ADMno1::gasPhaseRate()
{

    GRPtrs_[0] = para_.DTOS() * para_.kLa() 
               * (YPtrs_[7].internalField() - R_ * TopDummy_.internalField() * GPtrs_[0].internalField() * KHh2_);

    GRPtrs_[1] = para_.DTOS() * para_.kLa() 
               * (YPtrs_[8].internalField() - R_ * TopDummy_.internalField() * GPtrs_[1].internalField() * KHch4_);

    GRPtrs_[2] = para_.DTOS() * para_.kLa() // Sco2 instead of SIC
               * (MPtrs_[0].internalField() - R_ * TopDummy_.internalField() * GPtrs_[2].internalField() * KHco2_);
}


void Foam::ADMno1::gasSourceRate()
{

    // field of cell volume for mesh 
    scalarField volMeshField = GPtrs_[0].mesh().V().field();            

    // particle scaled gas volume
    scalarField volGas = volMeshField / (1.0 + (1.0/Vfrac_));

    // particle scaled liquid volume
    scalarField volLiq = volMeshField / (1.0 + Vfrac_);

    // particle scaled pipe resistance
    volScalarField kp
    (
        IOobject
        (
            "kp",
            fac_.mesh().time().timeName(),
            fac_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fac_.mesh(),
        dimensionedScalar
        (
           "kp_Default", 
            dimless,
            Zero
        )
    );

    // TODO: might lead to error if Gas dimension is different
    kp.dimensions().reset(dimVolume/dimTime/dimPressure);

    // TODO: actual volume would have effect on normalized kp
    kp.field() = para_.DTOS() * KP_ * (volMeshField / (Vgas_ + Vliq_).value()); // <---- this would be calculated

    //  volScalarField qGasLocal = kp * (Pgas - Pext_) * (Pgas / Pext_); 
    volScalarField qGasLocal = kp * (Pgas_ - Pext_);
    forAll( qGasLocal.field(), i )
    {
        if ( qGasLocal.field()[i] < 0.0 ) { qGasLocal.field()[i] = 1e-16; }
    }

    dGPtrs_[0].field() = GRPtrs_[0].field() * volLiq / volGas - GPtrs_[0].field() * qGasLocal.field() / volGas;

    dGPtrs_[1].field() = GRPtrs_[1].field() * volLiq / volGas - GPtrs_[1].field() * qGasLocal.field() / volGas;

    dGPtrs_[2].field() = GRPtrs_[2].field() * volLiq / volGas - GPtrs_[2].field() * qGasLocal.field() / volGas;
};


// ************************************************************************* //