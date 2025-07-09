[![Build status](https://github.com/WIAS-PDELib/ExtendableFEM.jl/workflows/linux-macos-windows/badge.svg)](https://github.com/WIAS-PDELib/ExtendableFEM.jl/actions)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://wias-pdelib.github.io/ExtendableFEM.jl/stable/index.html)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://wias-pdelib.github.io/ExtendableFEM.jl/dev/index.html)
[![DOI](https://zenodo.org/badge/668345991.svg)](https://zenodo.org/doi/10.5281/zenodo.10563834)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)

**ExtendableFEM.jl** is a high-level, extensible finite element method (FEM) library for Julia, supporting flexible problem descriptions, custom operators, and advanced grid management.

## Features

- High-level, extensible API for solving PDE problems by finite element methods
- Flexible `ProblemDescription` interface for assigning unknowns and operators
- Supports custom kernel functions for bilinear, linear, and nonlinear forms
- Automatic assembly and Newton's method for NonlinearOperators
- Builds upon [ExtendableGrids.jl](https://github.com/WIAS-PDELib/ExtendableGrids.jl) and low level structures from [ExtendableFEMBase.jl](https://github.com/WIAS-PDELib/ExtendableFEMBase.jl)

## Quick Example

The following minimal example demonstrates how to set up and solve a Poisson problem:

```julia
using ExtendableFEM
using ExtendableGrids

# Build a uniform-refined 2D unit square grid with triangles
xgrid = uniform_refine(grid_unitsquare(Triangle2D), 4)

# Create a new PDE description
PD = ProblemDescription()

# Define and assign the unknown
u = Unknown("u"; name = "potential")
assign_unknown!(PD, u)

# Assign Laplace operator (diffusion term)
assign_operator!(PD, BilinearOperator([grad(u)]; factor = 1e-3))

# Assign right-hand side data
function f!(fval, qpinfo)
    x = qpinfo.x # global coordinates of quadrature point
    fval[1] = x[1] * x[2]
end
assign_operator!(PD, LinearOperator(f!, [id(u)]))

# Assign Dirichlet boundary data (u = 0)
assign_operator!(PD, HomogeneousBoundaryData(u; regions = 1:4))

# Discretize: choose finite element space
FEType = H1Pk{1,2,3} # cubic H1-conforming element with 1 component in 2D
FES = FESpace{FEType}(xgrid)

# Solve the problem
sol = solve!(PD, [FES])

# Plot the solution
using PyPlot
plot(id(u), sol; Plotter = PyPlot)
```

## Citing

If you use ExtendableFEM.jl in your research, please cite [this Zenodo record](https://zenodo.org/doi/10.5281/zenodo.10563834).

## License
ExtendableFEM.jl is licensed under the MIT License. See [LICENSE](https://github.com/WIAS-PDELib/ExtendableFEM.jl/blob/master/LICENSE) for details.
