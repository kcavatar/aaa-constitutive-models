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

#include "incompressibleRaghavanVorpElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleRaghavanVorpElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, incompressibleRaghavanVorpElastic, nonLinGeomMechLaw
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::incompressibleRaghavanVorpElastic::incompressibleRaghavanVorpElastic
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
    mu_(2*alpha_),
    K_(dict.lookup("K"))
{
    if (planeStress())
    {
        notImplemented
        (
            type() + " mechanical law is not implemented for planeStress"
        );
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::incompressibleRaghavanVorpElastic::~incompressibleRaghavanVorpElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::incompressibleRaghavanVorpElastic::rho() const
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


Foam::tmp<Foam::volScalarField> Foam::incompressibleRaghavanVorpElastic::impK() const
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


void Foam::incompressibleRaghavanVorpElastic::correct
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

    // Calculate the hydrostatic stress
    updateSigmaHyd
    (
        0.5*K_*(pow(J, 2.0) - 1.0),
        (4.0/3.0)*mu_ + K_
    );

    // Calculate the left Cauchy-Green deformation tensor
    const volSymmTensorField B = symm(F() & F().T());

    // Compute invariants fields
    const volScalarField I1 = tr(B);
    // const volScalarField I2 = 0.5*(pow(tr(B), 2.0) - tr(B & B));

    // Calculate the deviatoric part of the Cauchy stress
    const volSymmTensorField s =
        2*
        (   alpha_
          + 2*beta_*(I1 - 3)
        )*B;

    // Calculate the Cauchy stress
    // for the Raghavan-Vorp model
    // The last RHS term is the initial hydrostatic pressure field
    // This term is important to assure the underformed configuration
    // to be stress-free
    sigma = sigmaHyd()*I + s - mu_*I;  // Initial hydrostatic pressure, p_0 = mu
}


void Foam::incompressibleRaghavanVorpElastic::correct(surfaceSymmTensorField& sigma)
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

    // Calculate pressure field with bulk modulus to approximate
    // incompressibility
    // Note: updateSigmaHyd is not used
    const surfaceScalarField sigmaHydf = 0.5*K_*(pow(J, 2.0) - 1.0);

    // Calculate the left Cauchy-Green deformation tensor
    const surfaceSymmTensorField Bf = symm(Ff() & Ff().T());

    // Compute invariants fields
    const surfaceScalarField I1 = tr(Bf);
    //const surfaceScalarField I2 = 0.5*(pow(tr(Bf), 2.0) - tr(Bf & Bf));

    // Calculate the deviatoric part of the Cauchy stress
    const surfaceSymmTensorField s =
        2*
        (   alpha_
          + 2*beta_*(I1 - 3)
        )*Bf;

    // Calculate the Cauchy stress
    // for the Mooney-Rivlin model
    // The last RHS term is the initial hydrostatic pressure field
    // This term is important to assure the underformed configuration
    // to be stress-free
    sigma = sigmaHydf*I + s - mu_*I; // Initial hydrostatic pressure, p_0 = mu
}


// ************************************************************************* //
