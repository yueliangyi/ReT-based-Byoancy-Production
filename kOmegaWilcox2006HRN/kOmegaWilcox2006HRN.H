/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Class
    Foam::RASModels::kOmegaWilcox2006HRN

Group
    grpRASTurbulence

Description
    Standard high Reynolds-number k-omega turbulence model for
    incompressible and compressible flows.

    References:
    \verbatim
        Wilcox, D. C. (1998).
        Turbulence modeling for CFD
        (Vol. 2, pp. 103-217). La Canada, CA: DCW industries.
    \endverbatim

    The default model coefficients are
    \verbatim
        kOmegaWilcox2006HRNCoeffs
        {
            Cmu         0.09;  // Equivalent to betaStar
            alpha       0.52;
            beta        0.072;
            alphak      0.5;
            alphaOmega  0.5;
        }
    \endverbatim

SourceFiles
    kOmegaWilcox2006HRN.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaWilcox2006HRN_H
#define kOmegaWilcox2006HRN_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "ReynoldsStress.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kOmegaWilcox2006HRN Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaWilcox2006HRN
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{

protected:

    // Protected data

            Switch transitionalBuoyancy_;

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar beta0_;
            dimensionedScalar gamma_;
            dimensionedScalar alphaK_;
            dimensionedScalar alphaOmega_;
            dimensionedScalar sigmad0_;

            dimensionedScalar alphaBInf_;
            dimensionedScalar Rbeta_;
            dimensionedScalar Rk_;

            dimensionedScalar Clim1_;
            dimensionedScalar Clim2_;


        // Fields

            volScalarField GbyNu_;
            volScalarField BbyNu_;
            volScalarField BK_;
            volScalarField k_;
            volScalarField omega_;
            volScalarField omegaTilde1_;
            volScalarField omegaTilde2_;
            uniformDimensionedVectorField g_;

            volScalarField pFlowIndex_;

            volScalarField beta_;
            volScalarField ReT_;
            volScalarField alphaB_;
            volScalarField alphaStar_;


    // Protected Member Functions

        virtual void correctNut(const volScalarField& S2);
        virtual void correctNut();
        virtual void correctB();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaWilcox2006HRN");


    // Constructors

        //- Construct from components
        kOmegaWilcox2006HRN
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kOmegaWilcox2006HRN() = default;


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //-
        virtual tmp<volScalarField> alphaStar() const
        {
            return
                (beta0_/3.0+ReT_/Rk_)
              / (1.0+ReT_/Rk_);
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    alphaK_*k_/omega_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DomegaEff",
                    alphaOmega_*k_/omega_ + this->nu()
                )
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the modified turbulence specific dissipation rate
        virtual tmp<volScalarField> omegaTilde1(const volScalarField& S2);
        virtual tmp<volScalarField> omegaTilde2();

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    Cmu_*k_*omega_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Return the turbulence strength
        virtual tmp<volScalarField> ReT() const
        {
            return k_/omega_/this->nu();
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmegaWilcox2006HRN.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
