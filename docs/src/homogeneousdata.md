# HomogeneousData

`HomogeneousData` provides a convenient way to enforce homogeneous (zero) Dirichlet boundary conditions or constraints in a finite element problem. It automatically sets the solution to zero on specified regions or entities.

## API Reference

```@autodocs
Modules = [ExtendableFEM]
Pages = ["common_operators/homogeneousdata_operator.jl"]
Order   = [:type, :function]
```

## Example: Zero Dirichlet Boundary Condition

```julia
# Impose u = 0 on boundary region 1
assign_operator!(PD, HomogeneousBoundaryData(u; regions = [1]))
```
