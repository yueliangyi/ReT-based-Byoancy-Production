# $Re_T$-based Byoancy Production

The buoyancy production term has been recognized as essential in turbulence closure models to address the persistent overestimation of turbulence near the air–water interface in Reynolds-averaged Navier–Stokes (RANS) simulations of breaking waves. While generally effective, two-equation k–ω based turbulence models typically use a simple representation of the buoyancy-production term with a constant closure coefficient in the turbulent kinetic energy equation. This can cause turbulence levels to collapse near the interface, essentially to zero within numerical accuracy, resulting in a loss of coupling between the two turbulence model equations. Such a breakdown inhibits the accurate initiation and evolution of turbulence, particularly during wave breaking onset. In this study, we address the decoupling problem by introducing a variable, turbulence-Reynolds-number-based closure coefficient for the buoyancy-production term. This adaptive formulation directly relates the strength of buoyancy production to local turbulence levels. For simulating spilling breaking waves, period-averaged surface elevation profiles show better agreement with experimental measurements. Wave-to-wave variability analysis further highlights the stabilizing effect of the proposed formulation. The notably improved undertow predictions support the use of this variable-coefficient approach in future RANS simulations of surface waves.

# Reference

Yue, L., and Li, Y., 2026. On the buoyancy production term for Reynolds-averaged modelling of breaking waves. Coastal Engineering 205, 104935.
DOI: https://doi.org/10.1016/j.coastaleng.2025.104935

# Installation

Download the repository
```
git clone https://github.com/yueliangyi/ReT-based-Byoancy-Production
```

Create folder for turbulence model (if the folders already exist skip this part)
```
mkdir -p $WM_PROJECT_USER_DIR/src/
```

Move the folder to the user source code
```
mv TurbulenceModels $WM_PROJECT_USER_DIR/src/
```

Go to the directory and compile the turbulence models
```
cd $WM_PROJECT_USER_DIR/src/TurbulenceModels
wmake libso
```



