# Parallel Assembly

Efficient assembly of finite element operators is crucial for large-scale simulations. This package supports parallel assembly of most operators.

## How It Works

- When `parallel = true` is set for an operator, assembly is performed in parallel over grid partitions (color groups).
- The computational grid must provide partition information. See the [ExtendableGrids.jl documentation on partitioning](https://wias-pdelib.github.io/ExtendableGrids.jl/stable/partitioning/) for details on how to partition your grid.
- The system matrix is constructed using a parallel-aware type (`MTExtendableSparseMatrixCSC` from [ExtendableSparse.jl](https://github.com/WIAS-PDELib/ExtendableSparse.jl)), which allows safe concurrent writes from different partitions.
- The solver automatically detects and uses the appropriate matrix type when the grid is partitioned.

## Usage

- Set the `parallel = true` keyword argument in the operators or in the solver configuration.
- Ensure your grid is partitioned appropriately for your problem size and hardware.
- For DG operators that assemble along cell faces, use the option `edges = true` in the partition call to enable parallel face assembly. Otherwise, face-based assembly will remain sequential.

## Example

```julia
# Partition the grid for parallel assembly
partition!(xgrid; nparts = 12, edges = true)

# Construct operator with parallel assembly enabled
assign_operator!(PD, BilinearOperator([grad(u)]; parallel = true))
```
