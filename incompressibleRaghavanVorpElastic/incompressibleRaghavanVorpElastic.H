/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    incompressibleRaghavanVorpElastic

Description
    Incompressible 2-parameter Raghavan-Vorp model, with strain-energy
    function defined as:

        psi = alpha*(I1 - 3) + beta*(I1 - 3)^2

    where I1 is the first invariant of the right Cauchy-Green deformation 
    tensor; alpha and beta are material constants.

    The Cauchy stress tensor is given by:

        sigma = -p*I + 2*[alpha + 2*beta*(I1 - 3)]*B
              = -p*I + 2*alpha*[1 + 2*(beta/alpha)*(I1 - 3)]*B
    =>  sigma = -p*I + mu*[1 + gamma*(I1 - 3)]*B

    where
    p        hydrostatic pressure (sigmaHyd == -p)
    I        identity tensor
    B	     left Cauchy-Green deformation tensor (== F & F.T)
    mu       2*alpha
    gamma    2*(beta/alpha)

    Incompressibility is enforced using the penalty method, where the bulk
    modulus (penalty) term should be set orders of magnitude greater than the
    shear modulus parameters.

    The bulk modulus acts as a penalty parameter, assuring that the
    incompressibility condition is satisfied. To ensure this, its value
    should be 'large enough' relative to the shear modulus parameters, but how
    large is subject to numerical experimentation. References below recommend it
    to be around 1e3 to 1e7 times the largest material parameter of the model
    (alpha or beta). Based on numerical experiments, a value of 1e3 times the
    small-strain shear modulus of the Mooney-Rivlin model (mu = 2*(c01 + c10))
    should be enough for the current implementation.

    Optionally: a Poisson's equation can be solved to calculate sigmaHyd. This
    is required as incompressibility is approached.
    (Iago Lessa de Oliveira, Philip Cardiff, UCD.)

 
    References
    - Raghavan et al. (2000). Toward a biomechanical tool to evaluate rupture 
      potential of abdominal aortic aneurysm: identification of a finite strain
      constitutive model and evaluation of its applicability.
      Journal of Biomechanics, 33, 475-482.
      
    - Oliveira et al. (2020). Implementation and numerical verification of an 
      incompressible three-parameter Mooney-Rivlin model for large deformation
      of soft tissues. 15th OpenFOAM Workshop (OFW15), June 22-25, 2020,
      Washington DC, USA
     
     
    - Bijelonja et al. (2005). A ﬁnite volume method for large strain analysis
      of incompressible hyperelastic materials. Int. J. Numer. Meth. Engg 2005;
      64:1594–1609 https://onlinelibrary.wiley.com/doi/epdf/10.1002/nme.1413
     
    - Nonlinear Solid Mechanics - A continuum approach for engineering
      Gerhard Hozapfel.


SourceFiles
    incompressibleRaghavanVorpElastic.C

Author
    G R Krishna Chand Avatar. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef raghavanVorpElastic_H
#define raghavanVorpElastic_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class incompressibleRaghavanVorpElastic
:
    public mechanicalLaw
{
    // Private data

        // Initial density
        const dimensionedScalar rho_;

        // First material parameter
        const dimensionedScalar alpha_;

        // Second material parameter
        const dimensionedScalar beta_;

        // Shear modulus: mu = 2*alpha
        const dimensionedScalar mu_; 

        // Bulk modulus
        const dimensionedScalar K_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        incompressibleRaghavanVorpElastic(const incompressibleRaghavanVorpElastic&);

        //- Disallow default bitwise assignment
        void operator=(const incompressibleRaghavanVorpElastic&);

public:

    //- Runtime type information
    TypeName("incompressibleRaghavanVorpElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        incompressibleRaghavanVorpElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~incompressibleRaghavanVorpElastic();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
