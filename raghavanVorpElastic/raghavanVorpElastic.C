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

\*---------------------------------------------------------------------------*/

#include "raghavanVorpElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(raghavanVorpElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, raghavanVorpElastic, nonLinGeomMechLaw
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::raghavanVorpElastic::raghavanVorpElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    rho_(dict.lookup("rho")),
    alpha_(dict.lookup("alpha")),
    beta_(dict.lookup("beta")),
    nu(dict.lookup("nu")),
    mu_(2*alpha_),
    K_("K", dimPressure, 0.0)
{
    
    // Define 2nd Lame parameter
    const dimensionedScalar lambda_ = (2.0*nu*mu_)/(1.0 - 2.0*nu);
    K_ = (2.0/3.0)*mu_ + lambda_;
    
    // Warning if plane stress analysis is chosen
    if (planeStress())
    {
        notImplemented
        (
            type() + " mechanical law is not implemented for planeStress"
        );
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::raghavanVorpElastic::~raghavanVorpElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::raghavanVorpElastic::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rhoLaw",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            calculatedFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::raghavanVorpElastic::impK() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            (4.0/3.0)*mu_ + K_
        )
    );
}


void Foam::raghavanVorpElastic::correct
(
    volSymmTensorField& sigma
)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J = det(F());

    // Calculate the volume preserving left Cauchy Green deformation tensor
    const volSymmTensorField B_iso = pow(J, -2.0/3.0)*symm(F() & F().T());

    // Compute invariant field I1 
    const volScalarField I1 = tr(B_iso);
    // const volScalarField I2 = 0.5*(pow(tr(B), 2.0) - tr(B & B));

    // Calculate the deviatoric part of the Cauchy stress tensor
    const volSymmTensorField s = dev(2.0*(alpha_ + 2.0*beta_*(I1 - 3))*B_iso);

    // Calculate the Cauchy stress for Raghavan-Vorp model
    sigma = (1.0/J)*(0.5*K_*(pow(J, 2) - 1)*I + s);
}


void Foam::raghavanVorpElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateFf(sigma, mu_, K_))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J = det(Ff());

    // Calculate the volume preserving left Cauchy Green deformation tensor
    const surfaceSymmTensorField B_iso = pow(J, -2.0/3.0)*symm(Ff() & Ff().T());

    // Compute invariants fields
    const surfaceScalarField I1 = tr(B_iso);
    //const surfaceScalarField I2 = 0.5*(pow(tr(Bf), 2.0) - tr(Bf & Bf));

    // Calculate the deviatoric part of the Cauchy stress
    const surfaceSymmTensorField s = dev(2.0*(alpha_ + 2.0*beta_*(I1 - 3))*B_iso);
    
    // Calculate the Cauchy stress for Raghavan-Vorp model
    sigma = (1.0/J)*(0.5*K_*(pow(J, 2) - 1)*I + s);
}


// ************************************************************************* //
