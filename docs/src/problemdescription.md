# Problem Description

The `ProblemDescription` is the central object for defining finite element problems in a high-level, flexible, and modular way. It encodes the weak form of your PDE, including unknowns, operators, and boundary or interface conditions, usually without requiring discretization details at this stage, but region numbers (for boundary conditions, etc.) may be referenced.


```@docs
ProblemDescription
```

## Constructors and Assignment Functions

Use the following functions to construct and modify a `ProblemDescription`:

```@autodocs
Modules = [ExtendableFEM]
Pages = ["problemdescription.jl"]
Order   = [:function]
```

## Unknowns

An `Unknown` represents a physical quantity (e.g., velocity, pressure, temperature) in the `ProblemDescription`. Unknowns are used to tag solution components and to specify which variables operators act on.

```@autodocs
Modules = [ExtendableFEM]
Pages = ["unknowns.jl"]
Order   = [:type, :function]
```

## Operators

Operators define how terms contribute to the system matrix or right-hand side. They can represent weak forms of differential operators, stabilization terms, constraints, or boundary conditions.

### Types of Operators

The main operator classes are:
- [NonlinearOperator](@ref) (e.g., nonlinear convection in Navier–Stokes)
- [BilinearOperator](@ref) (e.g., Laplacian in Poisson)
- [LinearOperator](@ref) (e.g., right-hand side in Poisson or Navier–Stokes)

For boundary conditions or global constraints, use:
- [InterpolateBoundaryData](@ref)
- [HomogeneousData](@ref)
- [FixDofs](@ref)
- [CombineDofs](@ref)

### Entities and Regions

Each operator assembles on certain mesh entities. The default is cell-wise assembly, but this can be changed via the `entities` keyword. Restrict assembly to subsets of entities using the `regions` keyword.

| Entities         | Description                                                      |
| :--------------- | :--------------------------------------------------------------- |
| AT_NODES         | Interpolate at mesh vertices (only for H1-conforming FEM)        |
| ON_CELLS         | Assemble/interpolate on mesh cells                               |
| ON_FACES         | Assemble/interpolate on all mesh faces                           |
| ON_IFACES        | Assemble/interpolate on interior mesh faces                      |
| ON_BFACES        | Assemble/interpolate on boundary mesh faces                      |
| ON_EDGES (*)     | Assemble/interpolate on all mesh edges   |
| ON_BEDGES (*)    | Assemble/interpolate on boundary mesh edges |

!!! note
    (*) = Only reasonable in 3D and still experimental; may have some issues.

### Function Operators

Function operators specify the mathematical operation to be applied to an unknown within an operator term. Each operator is defined as a pair of an `Unknown` (or integer index) and a `FunctionOperator` (such as `Identity`, `Gradient`, etc.), e.g., `(u, Identity)` or `(u, Gradient)`.

See the [full list of function operators](https://wias-pdelib.github.io/ExtendableFEMBase.jl/dev/functionoperators/).

For convenience and readability, common operator pairs have short aliases:
- `id(u)` for `(u, Identity)`
- `grad(u)` for `(u, Gradient)`
- `div(u)` for `(u, Divergence)`
- `curl1(u)` for `(u, CurlScalar)`
- `curl2(u)` for `(u, Curl2D)`
- `curl3(u)` for `(u, Curl3D)`
- `laplace(u)` or `Δ(u)` for `(u, Laplacian)`
- `hessian(u)` for `(u, Hessian)`
- `normalflux(u)` for `(u, NormalFlux)`
- `tangentialflux(u)` for `(u, TangentFlux)`
- `tangentialgrad(u)` for `(u, TangentialGradient)`
- `symgrad_voigt(u, factor)` for `(u, SymmetricGradient{factor})` in Voigt notation
- `εV(u, factor)` as a unicode alias for `symgrad_voigt(u, factor)`
- `grid(u)` for `(u, "grid")` (triggers gridplot in plot)
- `dofgrid(u)` for `(u, "dofgrid")` (triggers gridplot for the dofgrid in plot)
- `apply(u, FO)` for `(u, FO)` for any function operator type

Additional function operators for evaluating discontinuous quantities on faces (such as `jump`, `average`, `this`, `other`) are available.

For convenience, these discontinuous operator pairs also have short aliases:
- `jump((u, FO))` for `(u, Jump{FO})` (jump across a face)
- `average((u, FO))` for `(u, Average{FO})` (average across a face)
- `this((u, FO))` for `(u, Left{FO})` (value on the left side of a face)
- `other((u, FO))` for `(u, Right{FO})` (value on the right side of a face)

Here, `FO` can be any standard function operator (e.g., `Identity`, `Gradient`, etc.).
Of course, also something like `jump(id(u))` or `jump(grad(u))` works.
