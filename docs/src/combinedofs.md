# CombineDofs

`CombineDofs` provides a mechanism to couple degrees of freedom (DOFs) in a finite element system. This is especially useful for enforcing periodic boundary conditions, multi-point constraints, or other situations where DOFs from different locations or unknowns should be treated as identical in the assembled system.

```@autodocs
Modules = [ExtendableFEM]
Pages = ["common_operators/combinedofs.jl"]
Order   = [:type, :function]
```

## Periodic Boundary Conditions

To set up periodic boundary conditions, you often need to determine which DOFs should be coupled. The following utility functions can help:

```@docs
get_periodic_coupling_info
get_periodic_coupling_matrix
```

## Example: Periodic Boundary Coupling (extract from Example212)

Suppose you want to enforce periodicity between the left and right boundaries of a 2D domain. You can use `get_periodic_coupling_matrix` to find the corresponding DOFs, and then use `CombineDofs` to couple them in the system. The following code is adapted from [Example212_PeriodicBoundary2D.jl](https://github.com/WIAS-PDELib/ExtendableFEM.jl/blob/main/examples/Example212_PeriodicBoundary2D.jl):

```julia
# Define a function to map points on the left to the right boundary
function give_opposite!(y, x)
    y .= x
    y[1] = width - x[1]
    return nothing
end

# Compute the coupling matrix
coupling_matrix = get_periodic_coupling_matrix(FES, reg_left, reg_right, give_opposite!)

# Assign the CombineDofs operator to the problem description
assign_operator!(PD, CombineDofs(u, u, coupling_matrix))
```

See also [Example212](https://wias-pdelib.github.io/ExtendableFEM.jl/stable/examples/) and [Example312](https://wias-pdelib.github.io/ExtendableFEM.jl/stable/examples/) for more advanced use cases and details on periodic boundary conditions.
