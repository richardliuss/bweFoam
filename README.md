# bweFoam

OpenFOAM-based solver for depth intergrated Boussinesq-type equations.

This master branch has been tested OpenFOAM-7

When solving a problem, we recommend to use WENO-EXT(https://github.com/WENO-OF/WENOEXT) as the spacial derivates discretize scheme.


## Sofrware description:

### Governing equations

Solves the depth-intergrated 2D Boussinesq-type equations:

$$\frac{\partial h}{\partial t} + \nabla \cdot (h\vec{u}) = 0$$

$$\frac{\partial h\vec{u}}{\partial t} + (\vec{u} \cdot \nabla)h\vec{u} = -\left|g\right| h \nabla (z_{\rm{t}} + h) + \vec{\tau_{\rm{b}}} + \vec{A} + \vec{B}$$

where:

$h$: Flow depth

$\vec{u}$: Depth averaged velocity

$z_t$: Bottom elevation

$\vec{\tau_{b}}=\frac{gn_{\rm{c}}^{2}}{h^{1/3}}\left|\vec{u}\right|\vec{u}$: Bottom stresses

$\vec{A}=\frac{A_{\rm{m}}}{h}\nabla\cdot(h\nabla\vec{u})$: Eddy viscosity model

$\vec{B}=(1+\beta)\frac{h}{2}\nabla\left[\nabla\cdot\left(h\frac{\partial\vec{u}}{\partial t}\right)\right] + \beta\frac{h}{2}g\nabla\left[\nabla\cdot(h\nabla (z_{\rm{t}} + h))\right]  - (1+\beta)\frac{h^2}{6}\nabla(\nabla\cdot\frac{\partial\vec{u}}{\partial t}) - \beta\frac{h^2}{6}g\nabla\left[\nabla\cdot(\nabla (z_{\rm{t}} + h))\right]$

### Numerical schemes

This solver is tested using WENO-EXT(https://github.com/WENO-OF/WENOEXT) library. Therefore, we suggest you to use this library. We have forked it to this repo.

### Install the solver

1. Download the solver

2. Install WENO-EXT(optional)

3. Change to src folder

4. Type 'wmake'
