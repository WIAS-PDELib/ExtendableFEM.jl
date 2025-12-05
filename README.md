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
using ExtendableGrids # for grid management

# Build a structured ExtendableGrids.grid with mesh width h
h = 0.1
xgrid = simplexgrid(0:h:1, 0:h:1)

# Create a new PDE description
PD = ProblemDescription("Demo problem")

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

# Assign Dirichlet boundary data (u = 0) on all four grid boundary regions
assign_operator!(PD, HomogeneousBoundaryData(u; regions = [1,2,3,4]))

# Discretize: choose finite element space
FEType = H1Pk(order = 3, edim = 2, ncomponents = 1)
FES = FESpace{FEType}(xgrid)

# Solve the problem + unicode plot into terminal
sol = solve(PD, [FES]; plot = true, timeroutputs = :hide)

# Plot the solution
import PythonPlot
plot([id(u), grid(u)], sol; Plotter = PythonPlot)
```

## Running examples from documentation

In the [documentation](https://wias-pdelib.github.io/ExtendableFEM.jl/stable/index.html) many more examples can be found.
Each of these examples is implemented as a module that needs to be included first. Afterwards the main file of the module
can be run.
Usually the main function has a Plotter argument that can be used to toggle a plot of the solution with an
(already installed) backend of your choice (e.g. `Plotter=PythonPlot`, `GLMakie`, `Plots` or others supported by [GridVisualize.jl](https://github.com/WIAS-PDELib/GridVisualize.jl)).
Some examples need several further dependencies. To ensure an environment were every dependency is installed, one can use the test
environment via the package [TestEnv](https://github.com/JuliaTesting/TestEnv.jl). The following script runs Example201:
```julia
# activate test environment
using TestEnv
TestEnv.activate()

# include example file and load module
include("examples/Example201_PoissonProblem.jl")

# run example with default arguments
Example201_PoissonProblem.main()

# run with Plotting backend (added manually to the environment)
using PythonPlot
Example201_PoissonProblem.main(; Plotter = PythonPlot)
```


## Citing

If you use ExtendableFEM.jl in your research, please cite [this Zenodo record](https://zenodo.org/doi/10.5281/zenodo.10563834).

## License
ExtendableFEM.jl is licensed under the MIT License. See [LICENSE](https://github.com/WIAS-PDELib/ExtendableFEM.jl/blob/master/LICENSE) for details.
