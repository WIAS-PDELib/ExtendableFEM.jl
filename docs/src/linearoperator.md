# LinearOperator

`LinearOperator` provides a flexible interface for defining linear forms (right-hand side vectors) in finite element problems. These operators typically represent source terms, boundary data, or linearizations of nonlinear operators. The interface supports both standard and discontinuous Galerkin (DG) forms, and allows for custom assembly on cells, faces, or other grid entities.
It is also possible to assign a pre-computed vector as a LinearOperator.

## Constructors

```@autodocs
Modules = [ExtendableFEM]
Pages = ["common_operators/linear_operator.jl"]
Order   = [:type, :function]
```

## Example: Right-hand Side Operator

For a right-hand side operator of a Poisson problem with a given function `f(x)`, the kernel could look like:

```julia
function kernel!(result, qpinfo)
    result[1] = f(qpinfo.x)
end
u = Unknown("u")
LinearOperator(kernel!, [id(u)])
```

The argument `[id(u)]` specifies that the result of the kernel is multiplied with the identity evaluation of the test function.

## LinearOperatorDG

`LinearOperatorDG` is intended for linear forms that involve jumps or averages of discontinuous quantities on faces, requiring access to all degrees of freedom on neighboring cells.

```@autodocs
Modules = [ExtendableFEM]
Pages = ["common_operators/linear_operator_dg.jl"]
Order   = [:type, :function]
```
