# Stationary Solvers

This section describes how to solve stationary (time-independent) PDEs using the high-level API provided by this package. Both monolithic (single system) and block-iterative (subproblem) approaches are supported.
For time-dependent problems, there is an extra section.

## Meshes and FESpaces

To solve a `ProblemDescription`, you must provide discretization information:

- **Mesh:** The computational domain, represented as an `ExtendableGrid`. See [ExtendableGrids.jl](https://github.com/WIAS-PDELib/ExtendableGrids.jl) for details and grid constructors. For simplex grids, see [SimplexGridFactory.jl](https://github.com/WIAS-PDELib/SimplexGridFactory.jl). Gmsh mesh import is also supported.
- **Finite Element Spaces:** For each unknown, specify a finite element space (FESpace) that defines the ansatz functions. Construct with:

    ```julia
    FESpace{FEType}(grid::ExtendableGrid)
    ```
    where `FEType` is the finite element type. See the [list of available FETypes](https://wias-pdelib.github.io/ExtendableFEMBase.jl/dev/fems/) in the [ExtendableFEMBase.jl documentation](https://wias-pdelib.github.io/ExtendableFEMBase.jl/dev/).

## Monolithic Solve

To solve a problem, call `solve` with a `ProblemDescription` and an array of FESpaces (one per unknown). The solver automatically detects whether the problem is linear or nonlinear and chooses the appropriate algorithm (direct solve or fixed-point iteration/Newton's method). The nonlinear residual is reduced to the prescribed tolerance.

```@docs
solve
residual
```

!!! note
    The type of fixed-point algorithm depends on how nonlinearities are assembled. If all are assembled as `NonlinearOperator`, a Newton scheme is used (customizable via keyword arguments like `damping`). If nonlinearities are linearized by `LinearOperator` and `BilinearOperator`, other fixed-point iterations are used.

## Block-Iterative Solve (Subproblem Iteration)

For problems that can be solved by iterating over subproblems, configure each subproblem separately with a `ProblemDescription`/`SolverConfiguration`. This allows different tolerances and keyword arguments for each subproblem.

Construct a `SolverConfiguration` with:

```@docs
SolverConfiguration
```

Trigger the fixed-point iteration over subproblems with:

```@docs
iterate_until_stationarity
```
