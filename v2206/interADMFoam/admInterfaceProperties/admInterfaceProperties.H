
#ifndef admInterfaceProperties_H
#define admInterfaceProperties_H

#include "IOdictionary.H"
#include "surfaceTensionModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "interfaceProperties.H"

// namspace Foam?
namespace Foam
{

class admInterfaceProperties
:
    public interfaceProperties
{

// Private data
private:

    // Replicate of interfaceProperties.alpha1_
    const volScalarField& alpha1_;

    //- Sharp force coefficient (put 0.98-0.99 for static problems, 0.4-0.5 for dynamic)
    scalar cPc_;

    //-number of non-orthogonal corrector loop
    label nNonOrthogonalCorrectors_;
    
    // Capillary pressure
    volScalarField pc_;
    label pcRefCell_;
    scalar pcRefValue_;
    
    // Capillary flow    
    surfaceScalarField phic_;
        
    // Capillary pressure residual
    scalar eqnResidual_;

// Public data
public:

    // Constructors

    //- Construct from volume fraction field gamma and IOdictionary
    admInterfaceProperties
    (
        const volScalarField& alpha1,
        const volVectorField& U,
        const IOdictionary& dict
    );

    // Member functions 

    void calculatePhic();

    void solve();

    const volScalarField& pc()
    {
        return pc_;
    }

    const surfaceScalarField& phic()
    {
        return phic_;
    }
};

}

#endif