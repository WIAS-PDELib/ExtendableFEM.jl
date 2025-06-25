[![Build status](https://github.com/WIAS-PDELib/ExtendableFEM.jl/workflows/linux-macos-windows/badge.svg)](https://github.com/WIAS-PDELib/ExtendableFEM.jl/actions)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://wias-pdelib.github.io/ExtendableFEM.jl/stable/index.html)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://wias-pdelib.github.io/ExtendableFEM.jl/dev/index.html)
[![DOI](https://zenodo.org/badge/668345991.svg)](https://zenodo.org/doi/10.5281/zenodo.10563834)

# ExtendableFEM.jl

**ExtendableFEM.jl** is a high-level, extensible finite element method (FEM) library for Julia, designed for rapid prototyping and research in multiphysics and numerical analysis of FEM. It provides a flexible `ProblemDescription` interface, supports custom operators, and builds upon [ExtendableGrids.jl](https://github.com/WIAS-PDELib/ExtendableGrids.jl) and [ExtendableFEMBase.jl](https://github.com/WIAS-PDELib/ExtendableFEMBase.jl).

## Features

- High-level, modular API for PDE problems and finite element solvers
- Bilinear, linear, and nonlinear operators with custom kernel functions
- Automatic assembly and Newton's method for nonlinear operators
- Parallel assembly is possible


## Dependencies

- [ExtendableGrids.jl](https://github.com/WIAS-PDELib/ExtendableGrids.jl)
- [ExtendableSparse.jl](https://github.com/WIAS-PDELib/ExtendableSparse.jl)
- [ExtendableFEMBase.jl](https://github.com/WIAS-PDELib/ExtendableFEMBase.jl)
- [GridVisualize.jl](https://github.com/WIAS-PDELib/GridVisualize.jl)
- [DocStringExtensions.jl](https://github.com/JuliaDocs/DocStringExtensions.jl)
- [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)
- [DiffResults.jl](https://github.com/JuliaDiff/DiffResults.jl)

## Getting Started

The general workflow is as follows:

### 1. Geometry and Meshing

Create or load a grid using [ExtendableGrids.jl](https://github.com/WIAS-PDELib/ExtendableGrids.jl):

```julia
using ExtendableGrids
# Unit square, triangulated and refined
xgrid = uniform_refine(grid_unitsquare(Triangle2D), 4)
# Uniform rectangular grid
xgrid2 = simplexgrid(0:0.1:1, 0:0.2:2)
```

More complex grids can be created with [SimplexGridFactory.jl](https://github.com/WIAS-PDELib/SimplexGridFactory.jl) or by loading Gmsh files.

### 2. Problem Description

Pose your problem using the `ProblemDescription` interface:

```julia
PD = ProblemDescription()
u = Unknown("u"; name = "potential")
assign_unknown!(PD, u)
assign_operator!(PD, BilinearOperator([grad(u)]; factor = 1e-3))
f! = (result, qpinfo) -> (result[1] = qpinfo.x[1] * qpinfo.x[2])
assign_operator!(PD, LinearOperator(f!, [id(u)]))
assign_operator!(PD, HomogeneousBoundaryData(u; regions = 1:4))
```

### 3. Discretization

Choose finite element types and create spaces:

```julia
FES = FESpace{H1Pk{1,2,3}}(xgrid)
sol = FEVector(FES; tags = PD.unknowns)
fill(sol[u], 1) # optional: set initial values
```

### 4. Solve

Solve the problem:

```julia
sol = solve(PD, FES; init = sol)
```

For time-dependent problems, add time-derivative operators or use ODE solvers (see documentation).

### 5. Postprocessing and Plotting

After solving, postprocess or plot the solution:

```julia
using PyPlot
plot(id(u), sol; Plotter = PyPlot)
```

## Gradient-Robustness

ExtendableFEM.jl supports gradient-robust schemes via reconstruction operators or divergence-free elements. Gradient-robustness ensures that discrete velocities are independent of the exact pressure in incompressible flows. See the references below for more details.

## Citing

If you use ExtendableFEM.jl in your research, please cite [this Zenodo record](https://zenodo.org/doi/10.5281/zenodo.10563834).

## License

ExtendableFEM.jl is licensed under the MIT License. See [LICENSE](https://github.com/WIAS-PDELib/ExtendableFEM.jl/blob/main/LICENSE) for details.

## References

- [1] V. John, A. Linke, C. Merdon, M. Neilan, L. Rebholz, "On the divergence constraint in mixed finite element methods for incompressible flows", SIAM Review 59(3) (2017), 492--544. [Journal](https://doi.org/10.1137/15M1047696), [Preprint](http://www.wias-berlin.de/publications/wias-publ/run.jsp?template=abstract&type=Preprint&year=2015&number=2177)
- [2] A. Linke, C. Merdon, "Pressure-robustness and discrete Helmholtz projectors in mixed finite element methods for the incompressible Navier--Stokes equations", Computer Methods in Applied Mechanics and Engineering 311 (2016), 304--326. [Journal](http://dx.doi.org/10.1016/j.cma.2016.08.018), [Preprint](http://www.wias-berlin.de/publications/wias-publ/run.jsp?template=abstract&type=Preprint&year=2016&number=2250)
- [3] M. Akbas, T. Gallouet, A. Gassmann, A. Linke, C. Merdon, "A gradient-robust well-balanced scheme for the compressible isothermal Stokes problem", Computer Methods in Applied Mechanics and Engineering 367 (2020). [Journal](https://doi.org/10.1016/j.cma.2020.113069), [Preprint](https://arxiv.org/abs/1911.01295)
