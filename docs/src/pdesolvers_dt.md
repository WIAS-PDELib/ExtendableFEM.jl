# Time-dependent Solvers

This section describes how to solve time-dependent (non-stationary) PDEs using the high-level API.

## Approaches to Time-dependent Problems

- **Manual approach:**
    - Add custom time-derivative terms to the problem (e.g., a mass matrix as a `BilinearOperator` and `LinearOperator`s for previous time steps).
    - If more than one previous time step is needed (e.g., for BDF2 or multi-step methods), you must manage the storage and update of previous solutions manually.
- **Automatic approach:**
    - Reframe the `ProblemDescription` as an ODE problem and solve it using [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) via the `ExtendableFEMDiffEQExt.jl` extension.

Several time-dependent examples are available, including both approaches. See, for example, [Example103 (Burgers' equation)](https://wias-pdelib.github.io/ExtendableFEM.jl/stable/examples/) and [Example205 (Heat equation)](https://wias-pdelib.github.io/ExtendableFEM.jl/stable/examples/).

## Using SciMLBase.ODEProblem and DifferentialEquations.jl

You can reframe the spatial part of your PDE as the right-hand side of an `ODEProblem`. The `ProblemDescription` then describes the spatial operators and right-hand side:

```math
\begin{aligned}
M u_t(t) & = b(u(t)) - A(u(t)) u(t)
\end{aligned}
```

where:
- `A` and `b` are the assembled (linearized) spatial operator and right-hand side operators in the `ProblemDescription` (note: `A` comes with a minus sign).
- `M` is the mass matrix, which can be customized (as long as it remains constant).
- Operators in the `ProblemDescription` may depend on time (e.g., if their kernels use `qpinfo.time`) and will be reassembled at each time step.
- To avoid unnecessary reassembly, use `store = true` for operators that do not change in time, or `constant_matrix = true` in the `SolverConfiguration` to skip full matrix reassembly.

```@autodocs
Modules = [ExtendableFEM]
Pages = ["solvers_diffeq.jl", "diffeq_interface.jl"]
Order   = [:type, :function]
```

!!! note
    When using DifferentialEquations.jl, set `autodiff=false` in the solver options, as automatic differentiation of the generated ODEProblem with respect to time is not currently supported.

## Example: 2D Heat Equation (extracted from Example205)

The following `ProblemDescription` yields the space discretization of the heat equation (including homogeneous boundary conditions; equivalent to the Poisson equation):

```julia
PD = ProblemDescription("Heat Equation")
u = Unknown("u"; name = "temperature")
assign_unknown!(PD, u)
assign_operator!(PD, BilinearOperator([grad(u)]; store = true, kwargs...))
assign_operator!(PD, HomogeneousBoundaryData(u))
```

Given a finite element space `FES` and an initial `FEVector` `sol` for the unknown, the `ODEProblem` for a time interval `(0, T)` can be generated and solved as follows:

```julia
prob = generate_ODEProblem(PD, FES, (0, T); init = sol)
DifferentialEquations.solve(prob, ImplicitEuler(autodiff = false), dt = 1e-3, dtmin = 1e-6, adaptive = true)
```

## Tips

- For more advanced time-stepping schemes, manage previous solutions and time-derivative terms manually in the `ProblemDescription`.
- See the [examples](https://wias-pdelib.github.io/ExtendableFEM.jl/stable/examples/) for practical implementations of time-dependent problems.
- For further details on the ODE interface, see the [ExtendableFEMDiffEQExt.jl documentation](https://github.com/WIAS-PDELib/ExtendableFEMDiffEQExt.jl).
