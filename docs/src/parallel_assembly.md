# Parallel Assembly

Within the solve call, all operators will be assembled according to their configured parameters.
When 'parallel = true' is used (not available yet for all DG operators), the operators will be assembled in parallel
based on the color partitions within the grid. Hence, the computational grid must provide these
partitions, see [Documentation of ExtendableGrids.jl on Partitioning](https://wias-pdelib.github.io/ExtendableGrids.jl/stable/partitioning/)
for details. Also the sparse system matrix needs to be able to work on different
partitions in parallel. Once, the grid has partitions, the solver automatically
uses a suitable constructor for the system matrix (MTExtendableSparseMatrixCSC from [ExtendableSparse.jl](https://github.com/WIAS-PDELib/ExtendableSparse.jl)).


!!! note

    DG operators that assemble along cell faces need the option 'edges = true' in the partition call for the grid partitioning, otherwise assembly will be still sequentially.
    
