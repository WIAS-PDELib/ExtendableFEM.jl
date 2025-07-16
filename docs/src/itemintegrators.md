# Item Integrators

`ItemIntegrator` provides a flexible interface for computing derived quantities from finite element solutions. These include a posteriori error estimators, norms, physical quantities (e.g., drag/lift coefficients), and other statistics that are computed by integrating over mesh entities (cells, faces, etc.).

## API Reference

```@autodocs
Modules = [ExtendableFEM]
Pages = ["common_operators/item_integrator.jl"]
Order   = [:type, :function]
```

## ItemIntegratorDG

`ItemIntegratorDG` is intended for quantities that involve jumps or averages of discontinuous quantities on faces, requiring access to all degrees of freedom on neighboring cells. This is essential for DG methods and certain error estimators.

```@autodocs
Modules = [ExtendableFEM]
Pages = ["common_operators/item_integrator_dg.jl"]
Order   = [:type, :function]
```

See, e.g., Example207, Example210 or Example245 for some practical use cases.
