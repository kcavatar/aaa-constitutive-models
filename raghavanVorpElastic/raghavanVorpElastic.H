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
    raghavanVorpElastic

Description
    2-parameter Raghavan-Vorp (compressible) hyperelastic model for abdominal
    aortic aneurysms (AAA) with the decoupled form of its strain-energy 
    function is
        
        psi = psi_vol (J) + psi_iso (B_iso)
       
    where, the volumetric part is given by
           
        psi_vol (J) = 0.5*K*[0.5*(J^2 - 1) - log (J)]
         
    and the isochoric part is

        psi_iso (B_iso) = alpha*(I1 - 3) + beta*(I1 - 3)^2

    where
    
    alpha, beta   material constants
    J             Jacobian of the deformation gradient F
    B             left Cauchy–Green deformation tensor = F & F.T()
    B_iso         J^(-2.0/3.0)*B
    I1            first invariant of B_iso
    K             bulk modulus
    mu            shear modulus
    
    The Cauchy stress tensor is given by:

      sigma = 0.5*K*I*(J^2 - 1)/J + dev(2*[alpha + 2*beta*(I1 - 3)]*B_iso)/J
            = 0.5*K*I*(J^2 - 1)/J + dev(mu*[1 + (2*beta/alpha)*(I1 - 3)]*B_iso)/J

    where
    I        identity tensor
    B	     left Cauchy-Green deformation tensor (== F & F.T)
    mu       shear modulus (== 2.0*alpha)
    
    Relations:
    
    2nd Lame parameter, lambda = (2.0*nu*mu)/(1.0 - 2.0*nu); nu = Poisson's ratio
                            K  = (2.0/3.0)*mu + lambda

 
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
      Gerhard Hozapfel. (Page Nos: 229 - 233)




    @verbatim  Usage
    
    materialName
    {
        type            raghavanVorpElastic;
        rho             rho [1 -3 0 0 0 0 0]  1200;
        alpha           alpha   [1 -1 -2 0 0 0 0] 1.74e5; 
        beta            beta  [1 -1 -2 0 0 0 0]   1.881e6;  
        nu              nu  [0 0 0 0 0 0 0]   0.49; 
    }





SourceFiles
    raghavanVorpElastic.C

Author
    G R Krishna Chand Avatar. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef raghavanVorpElastic_H
#define raghavanVorpElastic_H

#include "mechanicalLaw.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "surfaceFields.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class raghavanVorpElastic
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
        
        // Poisson's ratio
        const dimensionedScalar nu;

        // Shear modulus: mu = 2*alpha
        dimensionedScalar mu_; 
        
        // Bulk modulus
        dimensionedScalar K_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        raghavanVorpElastic(const raghavanVorpElastic&);

        //- Disallow default bitwise assignment
        void operator=(const raghavanVorpElastic&);

public:

    //- Runtime type information
    TypeName("raghavanVorpElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        raghavanVorpElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~raghavanVorpElastic();


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
