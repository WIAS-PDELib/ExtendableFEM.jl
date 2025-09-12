# Restrictions

ExtendableFEM.jl provides functionality to apply various restrictions to your finite element problems through the `AbstractRestriction` type system. Restrictions are used to enforce additional constraints on your solution that cannot be easily expressed through standard boundary conditions or the weak formulation.

## Built-in Restrictions

### Mean Value Restriction

`MeanValueRestriction` allows you to constrain the mean value of a scalar unknown to a specific value. This is particularly useful in problems where solutions are only unique up to a constant.

```@autodocs
Modules = [ExtendableFEM]
Pages = ["common_restrictions/mean_value_restriction.jl"]
Order   = [:type, :function]
```

#### Example Usage
```julia
# Restrict the mean value of unknown u to 0
restriction = MeanValueRestriction(u)

# Restrict the mean value to 1.0 with a specific operator
restriction = MeanValueRestriction(u, value = 1.0, operator = MyCustomOperator)

# apply to the problem
assign_restriction!(PD, restriction)
```

### Coupled DOFs Restriction

`CoupledDofsRestriction` enables coupling between different degrees of freedom (DOFs) in your system. This is particularly useful for implementing periodic boundary conditions or other constraints where DOFs need to be related to each other. See also `CombineDofs` operator, which does the same by maniplulating the system matrix.
However, the application of `CoupledDofsRestriction` much faster for larger systems, since the manipulation of the system matrix, as performed by operators, is very expensive.

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
