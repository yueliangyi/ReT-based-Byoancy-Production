/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "kOmegaWilcox2006HRN.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaWilcox2006HRN<BasicTurbulenceModel>::correctNut
(
    const volScalarField& S2
)
{
    ReT_ = this->ReT();
    alphaStar_ = this->alphaStar();

    this->nut_ = k_/max(omegaTilde1(S2), omegaTilde2());

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void kOmegaWilcox2006HRN<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


template<class BasicTurbulenceModel>
void kOmegaWilcox2006HRN<BasicTurbulenceModel>::correctB()
{
    if (alphaBInf_.value() < SMALL) return;

    const volScalarField& rho =
        this->db().objectRegistry::lookupObject<volScalarField>("rho");
    volVectorField gradRhobyRho = fvc::grad(rho)/rho;

    alphaB_ = alphaBInf_;

    // Reference: Yue, L., & Li, Y. P. (2025)
    // On the buoyancy production term for Reynolds-averaged modelling of
    // breaking waves. Coastal Engineering, 104935.
    // DOI: 10.1016/j.coastaleng.2025.104935
    if (transitionalBuoyancy_)
    {
        alphaB_ *= pow(ReT_/Rbeta_,4)/(1.0 + pow(ReT_/Rbeta_,4));
    }

    BbyNu_ = alphaB_*(g_ & gradRhobyRho);

    //- Construct buoyancy production term for k
    // The buoyancy term should be affected by the stress limiters.
    BK_ = this->nut_*BbyNu_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaWilcox2006HRN<BasicTurbulenceModel>::kOmegaWilcox2006HRN
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    transitionalBuoyancy_
    (
        Switch::getOrAddToDict
        (
            "transitionalBuoyancy",
            this->coeffDict_,
            false
        )
    ),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    beta0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta0",
            this->coeffDict_,
            0.0708
        )
    ),
    gamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma",
            this->coeffDict_,
            0.52
        )
    ),
    alphaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            this->coeffDict_,
            0.6
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            this->coeffDict_,
            0.5
        )
    ),
    sigmad0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmad0",
            this->coeffDict_,
            0.125
        )
    ),

    alphaBInf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaBInf",
            this->coeffDict_,
            1.0/0.85
        )
    ),
    Rbeta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Rbeta",
            this->coeffDict_,
            8.0
        )
    ),
    Rk_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Rk",
            this->coeffDict_,
            6.0
        )
    ),

    Clim1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clim1",
            this->coeffDict_,
            0.0
        )
    ),
    Clim2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clim2",
            this->coeffDict_,
            0.05
        )
    ),

    GbyNu_
    (
        IOobject
        (
            IOobject::groupName("GbyNu", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Zero",dimless/pow(dimTime,2),Zero)
    ),
    BbyNu_
    (
        IOobject
        (
            IOobject::groupName("BbyNu", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Zero",dimless/pow(dimTime,2),Zero)
    ),
    BK_
    (
        IOobject
        (
            IOobject::groupName("BK", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Zero",dimViscosity/pow(dimTime,2),Zero)
    ),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omegaTilde1_
    (
        IOobject
        (
            IOobject::groupName("omegaTilde1", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Zero",dimless/dimTime,Zero)
    ),
    omegaTilde2_
    (
        IOobject
        (
            IOobject::groupName("omegaTilde2", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Zero",dimless/dimTime,Zero)
    ),
    g_
    (
        IOobject
        (
            "g",
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    pFlowIndex_
    (
        IOobject
        (
            IOobject::groupName("pFlowIndex", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0",dimless,0.0)
    ),

    beta_
    (
        IOobject
        (
            IOobject::groupName("beta", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0",dimless,0.0)
    ),

    ReT_
    (
        IOobject
        (
            IOobject::groupName("ReT", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Zero",dimless,Zero)
    ),
    alphaB_
    (
        IOobject
        (
            IOobject::groupName("alphaB", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Zero",dimless,Zero)
    ),
    alphaStar_
    (
        IOobject
        (
            IOobject::groupName("alphaStar", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Zero",dimless,Zero)
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);
    correctNut();

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaWilcox2006HRN<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        transitionalBuoyancy_.readIfPresent
        (
            "transitionalBuoyancy",
            this->coeffDict()
        );
        Cmu_.readIfPresent(this->coeffDict());
        beta0_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
        sigmad0_.readIfPresent(this->coeffDict());
        alphaBInf_.readIfPresent(this->coeffDict());
        Rbeta_.readIfPresent(this->coeffDict());
        Rk_.readIfPresent(this->coeffDict());
        Clim1_.readIfPresent(this->coeffDict());
        Clim2_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaWilcox2006HRN<BasicTurbulenceModel>::omegaTilde1
(
    const volScalarField& S2
)
{
    // Lambda1 limiter
    omegaTilde1_ = max
    (
        omega_,
        Clim1_*sqrt
        (
            max(S2, dimensionedScalar(S2.dimensions(),0.0))/Cmu_
        )
    );

    return omegaTilde1_;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaWilcox2006HRN<BasicTurbulenceModel>::omegaTilde2()
{
    tmp<volTensorField> tgradU = fvc::grad(this->U_);
    const volTensorField& gradU = tgradU();

    // Lambda2 limiter

    // Effectively potential ﬂow
    pFlowIndex_ = 2.0*magSqr(symm(gradU))/max
    (
        2.0*magSqr(skew(gradU)),
        dimensionedScalar("SMALL", dimless/(dimTime*dimTime), SMALL)
    );

    // Reference: Larsen BE, Fuhrman DR (2018)
    // On the over-production of turbulence beneath surface waves in
    // Reynolds-averaged Navier–Stokes models. Journal of Fluid Mechanics.
    // doi:10.1017/jfm.2018.577

    omegaTilde2_ = max
    (
        omega_,
        Clim2_*beta_/(Cmu_*gamma_)*pFlowIndex_*omega_
    );

    return omegaTilde2_;
}


template<class BasicTurbulenceModel>
void kOmegaWilcox2006HRN<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_) return;

    // Local references
    const alphaField& alpha = this->alpha_;
    const volVectorField& U = this->U_;
    const volScalarField& nut = this->nut_;
    const volScalarField& rho =
        this->db().objectRegistry::lookupObject<volScalarField>("rho");
    const surfaceScalarField& rhoPhi =
        this->db().objectRegistry::lookupObject<surfaceScalarField>("rhoPhi");

    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volTensorField& gradU = tgradU();
    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    // Compute coefficient beta
    volSymmTensorField S = symm(gradU);
    volSymmTensorField Shat = S - 1.0/2.0*tr(gradU)*I;
    volTensorField Omega = -skew(gradU);
    volScalarField absChi = mag(((Omega & Omega) && Shat)/pow(Cmu_*omega_,3));
    volScalarField fBeta = (1.0+85.0*absChi)/(1.0+100.0*absChi);
    beta_ = beta0_ * fBeta;

    // Compute coefficient sigmad
    volScalarField CDkOmega = fvc::grad(k_) & fvc::grad(omega_);
    volScalarField sigmad = pos(CDkOmega)*sigmad0_;


    GbyNu_ = fvc::grad(U) && dev(twoSymm(gradU));
    GbyNu_ = max
    (
        GbyNu_,
        dimensionedScalar("SMALL", GbyNu_.dimensions(), SMALL)
    );


    volScalarField G(this->GName(), nut*GbyNu_);

    // Compute buoyancy production terms
    correctB();

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();


    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alpha*rhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
        gamma_*alpha*rho*GbyNu_
      - fvm::SuSp(((2.0/3.0)*gamma_)*alpha()*rho()*divU(), omega_)
      - fvm::Sp(beta_()*alpha()*rho()*omega_(), omega_)
      - fvm::SuSp(-sigmad()/omega_()*alpha()*rho()*CDkOmega()/omega_(), omega_)
      + fvOptions(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alpha*rhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G()
      - fvm::SuSp(alpha()*rho()*BK_()/k_(), k_)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU(), k_)
      - fvm::Sp(Cmu_*alpha()*rho()*omega_(), k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(GbyNu_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
