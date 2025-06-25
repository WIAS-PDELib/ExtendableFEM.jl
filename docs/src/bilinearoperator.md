# BilinearOperator

`BilinearOperator` provides a flexible interface for defining bilinear forms (matrices) in finite element problems. These operators typically represent (linearizations of) PDE operators, stabilization terms, or custom user-defined couplings. The interface supports both standard and discontinuous Galerkin (DG) forms, and allows for custom assembly on cells, faces, or other grid entities.

## Overview

- **Standard Bilinear Operators:** Assemble contributions to the system matrix, including diffusion, mass, and convection terms, as well as custom couplings. Supports both cell- and face-based assembly, including jump terms for conforming and piecewise spaces that only require the dofs on that entity.
- **Discontinuous Galerkin (DG) Operators:** Use `BilinearOperatorDG` for forms that require evaluation of jumps or averages involving all degrees of freedom on neighboring cells (e.g., interior penalty stabilization, gradient jumps for H1 elements).
- **Custom Matrices:** You can also assign a user-assembled matrix as a `BilinearOperator`.

## Constructors

```@autodocs
Modules = [ExtendableFEM]
Pages = ["common_operators/bilinear_operator.jl"]
Order   = [:type, :function]
```

## BilinearOperatorDG

`BilinearOperatorDG` is intended for bilinear forms that involve jumps or averages of discontinuous quantities on faces, requiring access to all degrees of freedom on neighboring cells. This is essential for DG methods and certain stabilization techniques.

```@autodocs
Modules = [ExtendableFEM]
Pages = ["common_operators/bilinear_operator_dg.jl"]
Order   = [:type, :function]
```

## Examples

### Example: Stokes Operator

For the linear operator of a Stokes problem, a kernel could look like:

```julia
μ = 0.1 # viscosity parameter
function kernel!(result, input, qpinfo)
    ∇u, p = view(input,1:4), view(input, 5)
    result[1] = μ*∇u[1] - p[1]
    result[2] = μ*∇u[2]
    result[3] = μ*∇u[3]
    result[4] = μ*∇u[4] - p[1]
    result[5] = -(∇u[1] + ∇u[4])
    return nothing
end
```

The corresponding `BilinearOperator` constructor call reads:

```julia
u = Unknown("u"; name = "velocity")
p = Unknown("p"; name = "pressure")
BilinearOperator(kernel!, [grad(u), id(p)]; use_sparsity_pattern = true)
```

The `use_sparsity_pattern` argument ensures that the zero pressure-pressure block of the matrix is not assembled, since `input[5]` does not couple with `result[5]`.

### Example: Interior Penalty Stabilization

A popular stabilization for convection-dominated problems is based on jumps of the gradient, which can be realized with the following kernel:

```julia
function stab_kernel!(result, input, qpinfo)
    result .= input .* qpinfo.volume^2
end
```

and the `BilinearOperatorDG` constructor call:

```julia
u = Unknown("u")
assign_operator!(PD, BilinearOperatorDG(stab_kernel!, [jump(grad(u))]; entities = ON_IFACES, factor = 0.01))
```
