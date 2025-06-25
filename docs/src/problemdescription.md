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
- `NonlinearOperator` (e.g., nonlinear convection in Navier–Stokes)
- `BilinearOperator` (e.g., Laplacian in Poisson)
- `LinearOperator` (e.g., right-hand side in Poisson or Navier–Stokes)

For boundary conditions or global constraints, use:
- `InterpolateBoundaryData`
- `HomogeneousBoundaryData`
- `FixDofs`
- `CombineDofs`

### Entities and Regions

Each operator assembles on certain mesh entities. The default is cell-wise assembly, but this can be changed via the `entities` keyword. Restrict assembly to subsets of entities using the `regions` keyword.

| Entities         | Description                                                      |
| :--------------- | :--------------------------------------------------------------- |
| AT_NODES         | Interpolate at mesh vertices (only for H1-conforming FEM)        |
| ON_CELLS         | Assemble/interpolate on mesh cells                               |
| ON_FACES         | Assemble/interpolate on all mesh faces                           |
| ON_IFACES        | Assemble/interpolate on interior mesh faces                      |
| ON_BFACES        | Assemble/interpolate on boundary mesh faces                      |
| ON_EDGES (*)     | Assemble/interpolate on all mesh edges (3D only, experimental)   |
| ON_BEDGES (*)    | Assemble/interpolate on boundary mesh edges (3D only, experimental) |

!!! note
    (*) = Only reasonable in 3D and still experimental; may have some issues.

### Function Operators

Operators involve pairs of an `Unknown` and a `FunctionOperator` (or an alias such as `Identity`, `Gradient`, etc.). FunctionOperators specify the mathematical operation to be applied (see [the full list here](https://wias-pdelib.github.io/ExtendableFEMBase.jl/dev/functionoperators/)).

Additional FunctionOperators for evaluating discontinuous operators on faces are available (needed for DG methods or face terms in a posteriori error estimators):

```@autodocs
Modules = [ExtendableFEM]
Pages = ["jump_operators.jl"]
Order   = [:type, :function]
```
