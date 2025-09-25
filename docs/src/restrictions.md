# Restrictions

ExtendableFEM.jl provides functionality to apply various restrictions to your finite element problems through the `AbstractRestriction` type system. Restrictions enforce additional constraints on your solution that cannot be easily expressed through standard boundary conditions or the weak formulation.

## Built-in Restrictions

### Linear Functional Restriction

`LinearFunctionalRestriction` allows you to constrain a linear functional of a finite element unknown to a specific value. This is useful, for example, for enforcing that the solution has mean zero, or that an integral constraint is satisfied.

Documentation:
```@autodocs
Modules = [ExtendableFEM]
Pages = ["common_restrictions/linear_functional_restriction.jl"]
Order   = [:type, :function]
```

#### Example Usage
```julia
# Restrict the mean value of unknown u to 0
restriction = ZeroMeanValueRestriction(u)

# Restrict the mean value to 1.0 with a specific operator
restriction = MassRestriction(u, value = 1.0, operator = MyCustomOperator)

# Assign to the problem
assign_restriction!(PD, restriction)
```

### Coupled DOFs Restriction

`CoupledDofsRestriction` enables coupling between different degrees of freedom (DOFs) in your system. This is particularly useful for implementing periodic boundary conditions or other constraints where DOFs need to be related to each other. Compared to manipulating the system matrix directly via operators (e.g. `CombineDofs`), `CoupledDofsRestriction` is much faster for large systems.

Documentation:
```@autodocs
Modules = [ExtendableFEM]
Pages = ["common_restrictions/coupled_dofs_restriction.jl"]
Order   = [:type, :function]
```

#### Example Usage
```julia
# Create a coupling matrix for periodic boundary conditions
coupling_matrix = get_periodic_coupling_matrix(...)
restriction = CoupledDofsRestriction(coupling_matrix)
assign_restriction!(PD, restriction)
```

### Boundary Data Restriction

`BoundaryDataRestriction` enforces prescribed boundary values for a finite element unknown. It can be used for both homogeneous (zero) and non-homogeneous boundary conditions, interpolating the boundary data as needed.

Documentation:
```@autodocs
Modules = [ExtendableFEM]
Pages = ["common_restrictions/boundarydata_restriction.jl"]
Order   = [:type, :function]
```

#### Example Usage
```julia
# Homogeneous boundary data (default)
restriction = BoundaryDataRestriction(u)

# Inhomogeneous boundary data using a function
restriction = BoundaryDataRestriction(u, data = (result, qpinfo) -> result .= sin(qpinfo.x[1]))

assign_restriction!(PD, restriction)
```

## AbstractRestriction API

All restrictions are subtypes of `AbstractRestriction`. The following functions are available for interacting with restrictions:

- `assemble!(restriction, solution, SC; kwargs...)`: Assemble the internal matrices and vectors for the restriction.
- `restriction_matrix(restriction)`: Get the restriction matrix.
- `restriction_rhs(restriction)`: Get the right-hand side vector for the restriction.
- `fixed_dofs(restriction)`: Get the list of fixed degrees of freedom affected by the restriction.
- `is_timedependent(restriction)`: Returns whether the restriction is time-dependent.

---
