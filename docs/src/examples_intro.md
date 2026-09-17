# About the Examples

The examples in this package are designed to be practical, reproducible, and educational. They demonstrate a wide range of finite element applications and PDE model problems.

## Design Principles

- All examples can be run directly from the Julia REPL.
- Each example is a Julia module named after the file.
- Examples can serve as templates for your own projects.
- Many examples include test cases for automated verification.

## Running the Examples

To run an example (e.g., `Example212_PeriodicElasticity2D`):

1. Download the example file (see the source code link at the top of the example page).
2. Start Julia in the root directory of the package (or in the directory that contains
   the example file). Some examples need additional dependencies beyond the core
   requirements, so the easiest way to get all of them is to use the package
   [test environment](https://github.com/JuliaTesting/TestEnv.jl).
   After `Pkg.activate(".")`, it can be activated via the package
   [TestEnv](https://github.com/JuliaTesting/TestEnv.jl):

    ```julia
    julia> using Pkg; Pkg.activate(".")
    julia> using TestEnv; TestEnv.activate()
    ```

3. Include the example file and call its `main` function:

    ```julia
    julia> include("Example212_PeriodicElasticity2D.jl")
    julia> Example212_PeriodicElasticity2D.main()
    ```

4. Some examples offer visual output via the optional argument `Plotter = PyPlot` or `Plotter = GLMakie` (provided the package is installed and loaded):

    ```julia
    julia> Example212_PeriodicElasticity2D.main(Plotter = PyPlot)
    ```

## Overview of the Examples

The examples are listed below in numerical order. Each name links to the documentation
page of the example; the **keywords** indicate the application and the finite element
concepts demonstrated there (method type, coupling, error estimation, etc.).

| Example | Description |
|:---|:---|
| [Example103](module_examples/Example103_BurgersEquation.md) | 1D Burgers equation (scalar conservation law) with periodic boundary conditions, solved by a `DifferentialEquations`-based and a manual implicit Euler time discretization — **keywords:** time-dependent, nonlinear, periodic boundary |
| [Example105](module_examples/Example105_NonlinearPoissonEquation.md) | 1D nonlinear Poisson equation with Dirichlet data, solved by Newton iteration — **keywords:** conforming FEM, nonlinear |
| [Example106](module_examples/Example106_NonlinearDiffusion.md) | 1D porous-medium (nonlinear) diffusion equation with Neumann data — **keywords:** conforming FEM, time-dependent, nonlinear |
| [Example108](module_examples/Example108_RobinBoundaryCondition.md) | 1D convection–diffusion–reaction with a mixed Robin/Dirichlet boundary condition, verified against an exact solution — **keywords:** conforming FEM, Robin boundary condition, nonlinear, convergence study |
| [Example201](module_examples/Example201_PoissonProblem.md) | 2D Poisson problem on the unit square; demonstrates `BoundaryDataRestriction`- and operator-based boundary data — **keywords:** conforming FEM, restrictions |
| [Example202](module_examples/Example202_MixedPoissonProblem.md) | 2D Poisson problem in mixed (stress–potential) form with H(div)-conforming stress — **keywords:** conforming FEM, mixed method, H(div) |
| [Example203](module_examples/Example203_PoissonProblemDG.md) | 2D Poisson problem with inhomogeneous Dirichlet data, solved with the discontinuous Galerkin (interior penalty) method — **keywords:** discontinuous Galerkin, interior penalty |
| [Example204](module_examples/Example204_LaplaceEVProblem.md) | eigenvalue problem for the Laplacian on the L-shaped domain, using a KrylovKit iterative solver — **keywords:** conforming FEM, eigenvalue problem |
| [Example205](module_examples/Example205_HeatEquation.md) | 2D heat equation with Dirichlet data and time integration via `DifferentialEquations` — **keywords:** conforming FEM, time-dependent |
| [Example206](module_examples/Example206_CoupledSubGridProblems.md) | two Poisson problems on subdomains, coupled across a new internal interface region by a penalized interface condition — **keywords:** conforming FEM, coupled problem, interface conditions, penalty |
| [Example207](module_examples/Example207_AdvectionUpwindDG.md) | 2D advection equation with an upwind discontinuous Galerkin and inflow boundary data — **keywords:** discontinuous Galerkin, upwind stabilization, advection |
| [Example210](module_examples/Example210_LshapeAdaptivePoissonProblem.md) | residual error estimator and adaptive mesh refinement for the Poisson problem on the L-shaped domain; also demonstrates second-order derivative evaluation — **keywords:** conforming FEM, error estimation, adaptive mesh refinement, convergence study |
| [Example211](module_examples/Example211_LshapeAdaptiveEQPoissonProblem.md) | locally equilibrated (dual stress-reconstruction) error estimator, assembled in parallel on non-overlapping node patch groups — **keywords:** conforming + mixed FEM, equilibrated error estimation, adaptive mesh refinement, low-level assembly, convergence study |
| [Example212](module_examples/Example212_PeriodicElasticity2D.md) | linear elasticity on an unstructured periodic 2D grid with periodic coupling along one axis — **keywords:** conforming FEM, elasticity, periodic boundary, homogenization |
| [Example220](module_examples/Example220_ReactionConvectionDiffusion.md) | convection–diffusion–reaction; compares a conforming discretization with and without a gradient-jump stabilization for small diffusion — **keywords:** conforming FEM, interior penalty stabilization, convergence study |
| [Example225](module_examples/Example225_ObstacleProblem.md) | obstacle problem (Dirichlet energy minimization above an obstacle) via a penalty term and Newton iteration — **keywords:** conforming FEM, variational inequality, penalty, nonlinear |
| [Example226](module_examples/Example226_Thermoforming.md) | thermoforming benchmark: nonlinear elasticity with temperature-dependent material parameters — **keywords:** conforming FEM, nonlinear elasticity |
| [Example227](module_examples/Example227_ObstacleProblemLVPP.md) | the obstacle problem of Example225, solved with a latent-variable proximal-point algorithm instead of a penalty — **keywords:** conforming FEM, variational inequality, latent variable proximal point, nonlinear |
| [Example230](module_examples/Example230_NonlinearElasticity.md) | nonlinear (bimetal) elasticity; demonstrates region- and parameter-dependent nonlinear expressions — **keywords:** conforming FEM, nonlinear elasticity, convergence study |
| [Example235](module_examples/Example235_StokesIteratedPenalty.md) | Hagen–Poiseuille flow with the iterated penalty method (Bernardi–Raugel element) — **keywords:** incompressible flow, Stokes, iterated penalty method |
| [Example240](module_examples/Example240_SVRTEnrichment.md) | Stokes flow with a Scott–Vogelius type element: continuous P_k velocity functions enriched by Raviart–Thomas functions plus a discontinuous P_{k-1} pressure, giving a pointwise divergence-free velocity on general meshes, with DOF reduction — **keywords:** incompressible flow, Stokes, H(div) velocity enrichment, inf-sup stability, convergence study |
| [Example245](module_examples/Example245_NSEFlowAroundCylinder.md) | DFG benchmark (flow around a cylinder) on an externally generated mesh; computes drag/lift coefficients — **keywords:** incompressible flow, Navier–Stokes, nonlinear, external mesh generation, postprocessing |
| [Example250](module_examples/Example250_NSELidDrivenCavity.md) | Navier–Stokes in a lid-driven cavity over a cone; demonstrates vortex formation — **keywords:** incompressible flow, Navier–Stokes, nonlinear |
| [Example252](module_examples/Example252_NSEPlanarLatticeFlow.md) | planar lattice flow Navier–Stokes problem with exact data; compares the L2 error of two restriction-based discretizations — **keywords:** incompressible flow, Navier–Stokes, nonlinear, restrictions, convergence study |
| [Example260](module_examples/Example260_AxisymmetricNavierStokesProblem.md) | three-dimensional stagnation-point flow via the 2.5D axisymmetric Navier–Stokes formulation — **keywords:** incompressible flow, Navier–Stokes, axisymmetric |
| [Example264](module_examples/Example264_StokesDarcy.md) | coupled Stokes (free flow) / Darcy (porous medium) problem with interface transmission conditions — **keywords:** incompressible flow, Stokes + Darcy, interface conditions |
| [Example265](module_examples/Example265_FlowTransport.md) | Stokes flow in an Ω-shaped pipe driving a convection–diffusion transport equation; compares stabilized convection discretizations and uses a `CallbackOperator` — **keywords:** incompressible flow, Stokes, transport, upwind/grad-jump stabilization, callback operator |
| [Example270](module_examples/Example270_NaturalConvectionProblem.md) | natural convection (Boussinesq) in a triangular cavity, coupling Navier–Stokes and heat conduction — **keywords:** incompressible flow, Navier–Stokes, heat equation, convergence study |
| [Example275](module_examples/Example275_OptimalControlStokes.md) | optimal control for the Stokes problem; computes the optimal control with an adjoint formulation — **keywords:** incompressible flow, Stokes, optimal control, adjoint, convergence study |
| [Example280](module_examples/Example280_CompressibleStokes.md) | 2D compressible Stokes flow with an equation of state; a well-balanced (gradient-robust) discretization is compared against a standard one — **keywords:** incompressible flow, Stokes, compressible, equation of state, convergence study |
| [Example282](module_examples/Example282_IncompressibleMHD.md) | stationary incompressible viscous MHD: Navier–Stokes coupled with the magnetic induction equation — **keywords:** incompressible flow, Navier–Stokes, MHD, nonlinear |
| [Example284](module_examples/Example284_LevelSetMethod.md) | time-dependent convection of a level set by a given velocity field — **keywords:** conforming FEM, level set, free boundary, time-dependent |
| [Example285](module_examples/Example285_CahnHilliard.md) | phase-field Cahn–Hilliard equation in mixed form — **keywords:** mixed method, phase field, time-dependent |
| [Example290](module_examples/Example290_PoroElasticity.md) | three-field Biot consolidation model with an H(div)-conforming reconstruction to avoid Poisson locking — **keywords:** mixed method, poroelasticity, H(div), time-dependent, convergence study |
| [Example295](module_examples/Example295_SlidingDroplet.md) | droplet sliding down a surface with surface tension and slip, tracked with an ALE moving mesh up to a stationary state — **keywords:** incompressible flow, Navier–Stokes, ALE moving mesh, free boundary, time-dependent, nonlinear |
| [Example301](module_examples/Example301_PoissonProblem.md) | 3D Poisson problem on the unit cube; demonstrates iterative (Krylov) solvers with incomplete LU preconditioners — **keywords:** conforming FEM, 3D, Krylov solvers and preconditioning |
| [Example310](module_examples/Example310_DivFreeBasis.md) | best approximation of a divergence-free velocity by a divergence-free Raviart–Thomas (H(curl)-type) basis with a linearly independent basis construction — **keywords:** mixed method, H(curl), divergence-free basis, 3D |
| [Example312](module_examples/Example312_PeriodicElasticity3D.md) | 3D periodic elasticity on a generated unstructured grid with periodic coupling along one axis — **keywords:** conforming FEM, elasticity, periodic boundary, 3D |
| [Example313](module_examples/Example313_PeriodicPoisson.md) | 3D periodic Poisson problem verifying the restriction-based periodic boundary operator against the exactly linear solution — **keywords:** conforming FEM, periodic boundary, restrictions, 3D |
| [Example330](module_examples/Example330_HyperElasticity.md) | neo-Hookian hyperelasticity in 3D; the energy is twice differentiated automatically to set up the Newton scheme — **keywords:** conforming FEM, hyperelasticity, nonlinear, automatic differentiation, 3D |


## Learn a Concept

If you are looking for a particular finite element concept instead of a particular
problem, this table shows where to look: the corresponding documentation page(s) and
the examples that demonstrate the concept. Where no documentation page exists, the
listed examples are currently the primary reference.

| Concept | Documentation | Examples |
|:---|:---|:---|
| Nonlinear equations (Newton iteration) | [NonlinearOperator](nonlinearoperator.md) | [Example105](module_examples/Example105_NonlinearPoissonEquation.md), [Example225](module_examples/Example225_ObstacleProblem.md), [Example227](module_examples/Example227_ObstacleProblemLVPP.md), [Example230](module_examples/Example230_NonlinearElasticity.md), [Example330](module_examples/Example330_HyperElasticity.md) |
| Boundary conditions (Dirichlet, Neumann, Robin) | [Restrictions](restrictions.md), [Interpolate boundary data](interpolateboundarydata.md) | [Example108](module_examples/Example108_RobinBoundaryCondition.md), [Example201](module_examples/Example201_PoissonProblem.md), [Example313](module_examples/Example313_PeriodicPoisson.md) |
| Periodic boundary | [Restrictions](restrictions.md), [Combined DOFs](combinedofs.md) | [Example212](module_examples/Example212_PeriodicElasticity2D.md), [Example312](module_examples/Example312_PeriodicElasticity3D.md), [Example313](module_examples/Example313_PeriodicPoisson.md) |
| Variational inequalities (obstacle problem) | — | [Example225](module_examples/Example225_ObstacleProblem.md), [Example227](module_examples/Example227_ObstacleProblemLVPP.md) |
| Elasticity and hyperelasticity | — | [Example212](module_examples/Example212_PeriodicElasticity2D.md), [Example226](module_examples/Example226_Thermoforming.md), [Example230](module_examples/Example230_NonlinearElasticity.md), [Example312](module_examples/Example312_PeriodicElasticity3D.md), [Example330](module_examples/Example330_HyperElasticity.md) |
| Discontinuous Galerkin (nonconforming) | — | [Example203](module_examples/Example203_PoissonProblemDG.md), [Example207](module_examples/Example207_AdvectionUpwindDG.md) |
| Stabilization (interior penalty, upwind) | — | [Example203](module_examples/Example203_PoissonProblemDG.md), [Example207](module_examples/Example207_AdvectionUpwindDG.md), [Example220](module_examples/Example220_ReactionConvectionDiffusion.md), [Example265](module_examples/Example265_FlowTransport.md) |
| Mixed methods (H(div), H(curl)) | — | [Example202](module_examples/Example202_MixedPoissonProblem.md), [Example211](module_examples/Example211_LshapeAdaptiveEQPoissonProblem.md), [Example240](module_examples/Example240_SVRTEnrichment.md), [Example290](module_examples/Example290_PoroElasticity.md), [Example310](module_examples/Example310_DivFreeBasis.md) |
| Error estimation & adaptive refinement | — | [Example210](module_examples/Example210_LshapeAdaptivePoissonProblem.md), [Example211](module_examples/Example211_LshapeAdaptiveEQPoissonProblem.md) |
| Incompressible flow (Stokes, Navier–Stokes) | — | [Example235](module_examples/Example235_StokesIteratedPenalty.md), [Example240](module_examples/Example240_SVRTEnrichment.md), [Example245](module_examples/Example245_NSEFlowAroundCylinder.md), [Example250](module_examples/Example250_NSELidDrivenCavity.md), [Example252](module_examples/Example252_NSEPlanarLatticeFlow.md), [Example260](module_examples/Example260_AxisymmetricNavierStokesProblem.md), [Example264](module_examples/Example264_StokesDarcy.md), [Example265](module_examples/Example265_FlowTransport.md), [Example270](module_examples/Example270_NaturalConvectionProblem.md), [Example275](module_examples/Example275_OptimalControlStokes.md), [Example280](module_examples/Example280_CompressibleStokes.md), [Example282](module_examples/Example282_IncompressibleMHD.md), [Example295](module_examples/Example295_SlidingDroplet.md) |
| Quantities of interest (drag/lift, Nusselt number, etc.) | [Item integrators](itemintegrators.md) | [Example245](module_examples/Example245_NSEFlowAroundCylinder.md), [Example265](module_examples/Example265_FlowTransport.md), [Example270](module_examples/Example270_NaturalConvectionProblem.md), [Example295](module_examples/Example295_SlidingDroplet.md) |
| Coupled / multi-physics problems | — | [Example206](module_examples/Example206_CoupledSubGridProblems.md), [Example264](module_examples/Example264_StokesDarcy.md), [Example265](module_examples/Example265_FlowTransport.md), [Example270](module_examples/Example270_NaturalConvectionProblem.md), [Example282](module_examples/Example282_IncompressibleMHD.md), [Example290](module_examples/Example290_PoroElasticity.md) |
| Time-dependent problems | [Time-dependent solvers](pdesolvers_dt.md) | [Example103](module_examples/Example103_BurgersEquation.md), [Example106](module_examples/Example106_NonlinearDiffusion.md), [Example205](module_examples/Example205_HeatEquation.md), [Example284](module_examples/Example284_LevelSetMethod.md), [Example285](module_examples/Example285_CahnHilliard.md), [Example290](module_examples/Example290_PoroElasticity.md), [Example295](module_examples/Example295_SlidingDroplet.md) |
| Free boundary, level set, phase field, ALE | — | [Example284](module_examples/Example284_LevelSetMethod.md), [Example285](module_examples/Example285_CahnHilliard.md), [Example295](module_examples/Example295_SlidingDroplet.md) |
| Optimal control & adjoints | — | [Example275](module_examples/Example275_OptimalControlStokes.md) |
| Eigenvalue problems | — | [Example204](module_examples/Example204_LaplaceEVProblem.md) |
| Stationary solvers & preconditioning | [Stationary solvers](pdesolvers.md) | [Example301](module_examples/Example301_PoissonProblem.md) |
| Parallel assembly | [Parallel assembly](parallel_assembly.md) | [Example211](module_examples/Example211_LshapeAdaptiveEQPoissonProblem.md) |
| Postprocessing & visualization | [Postprocessing](postprocessing.md) | [Example210](module_examples/Example210_LshapeAdaptivePoissonProblem.md), [Example211](module_examples/Example211_LshapeAdaptiveEQPoissonProblem.md), [Example245](module_examples/Example245_NSEFlowAroundCylinder.md) |
