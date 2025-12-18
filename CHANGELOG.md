# CHANGES

## v1.8.0

### Added
  - new keyword streamlines that can be used in the plotfunction to trigger GridVisualize.streamplot, now used in Example252
  - new Example295 on a sliding droplet


## v1.7.2

### Fixed
  - fixed control flow and nonlinear residual computation in presence of restrictions

## v1.7.1

### Fixed
  - fixed type issues in LinearOperator and ItemIntegrator
  - fixed region awareness and output array sizing for boundary entities in ItemIntegratorDG

## v1.7.0

### Added
  - new parameter `autodiff_backend` for NonlinearOperator to change differentiation backend

### Changed
  - sparsity patterns of local jacobians are now resolved by SparseConnectivityTracer
  - local jacobians are now prepared by DifferentiationInterface


## v1.6.0 November 18, 2025

### Added

  - `CoupledDofsRestriction` accepts now an array of multiple coupling matrices.

### Fixed

  - Periodic coupling with multiple redundant coupling matrices are now reduced to full rank to make the system uniquely solvable.

## v1.5.0 September 15, 2025

### Added

  - `compute_periodic_coupling_matrix` can be computed thread parallel
  - new "Restriction" objects using Lagrange multipliers internally
  - user can apply linear `LinearFunctionalRestriction`, `ZeroMeanValueRestriction`, `MassRestriction`, `CoupledDofsRestriction` or `BoundaryDataRestriction` to the `ProblemDescription` by `apply_restriction!()`
  - Incorporated new restriction objects in examples 201, 252, 312 and 212.

## v1.4.0 July 23, 2025

### Changed
  - improved docstrings, documentation
  - improved print_convergencehistory and show functions

## v1.3.0 June 17, 2025

### Added
  - `plot` and `scalarplot!` for broken FE Spaces are now plotted discontinuously to avoid averaging of values at the grid nodes.
     This can be disabled by passing the kwarg `average_broken_plots = true` to the `plot` call.
  - new `id(u, c)` operator for the component `c` of a FE function `u` wrapping `IdentityComponent{c}` from ExtendableFEMBase

## v1.2.0 June 4, 2025

### Changed
  - TimerOutputs for measuring/storing/showing runtime and allocations in solve, now also for separate operators
  - Coupling matrix result of `compute_periodic_coupling_matrix` is no longer transposed
  - rewrote internals of `CombineDofs` `apply_penalty` method to speed up the assembly

### Fixed
  - HomogeneousData/InterpolateBoundaryData operator fix when system matrix is of type GenericMTExtendableSparseMatrixCSC

## v1.1.1 April 29, 2025

### Fixed
  - FixDofs operator does not crash when system matrix is of type GenericMTExtendableSparseMatrixCSC

## v1.1.0 April 17, 2025

### Changed
  - `get_periodic_coupling_matrix` with grid argument now deprecated

### Fixed
  - assembly times and allocations from solver times are now correctly measured

## v1.0.1 April 16, 2025

### Fixed

  - `assemble!` of (args-dependent) BilinearOperatorDG now works

## v1.0.0 April 14, 2025

### Added

  - new example on coupled Stokes-Darcy (Example264)

### Changed

  - example `Example301` now also demonstrates a nonlinear problem solved by an iterative linear solver
    with preconditioning.

  - `solve` uses now the residual equation for the linear systems
  - facelift `Example250`.

### Fixed

  - some bugfixes concerning finite elements on subregions and BilinearOperatorDG, related to issue #43

## v0.9.0 March 22, 2025

### Added

  - tests for explicit imports and undocumented names
  - Runic.jl code formatting
  - new operators `symgrad_Voigt` and `εV` for symmetric gradients in Voigt notation
  - new periodic boundary matrix assembler `get_periodic_coupling_matrix`
  - new `CombineDofs` operator accepting a periodic `coupling_matrix`

## v0.8 November 1, 2024
  - started changelog
  - first registered version since repository move (fixed some URL links in the doc and readme)
  - allow function-like kernel also in LinearOperator and BilinearOperator (thanks to @pjaap)

## October 28, 2024

Moved repository from https://github.com/chmerdon/ExtendableFEM.jl to https://github.com/WIAS-PDELib/ExtendableFEM.jl.
[WIAS-PDELib](https://github.com/WIAS-PDELib/) is a github organization created to collectively manage the Julia packages developed under
the lead of the [WIAS Numerical Mathematics and Scientific Computing](https://wias-berlin.de/research/rgs/fg3)  research group.
According to the [github docs on repository transfer](https://docs.github.com/en/repositories/creating-and-managing-repositories/transferring-a-repository#whats-transferred-with-a-repository),
all links to the previous repository location are automatically redirected to the new location, and all relationships with forks stay intact.
