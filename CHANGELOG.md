# CHANGES

## next version

## v1.0.0 April 11, 2025

### Added

  - new example on coupled Stokes-Darcy (Example264)

### Added

  - example `Example301` now also demonstrates a nonlinear problem solved by an iterative linear solver
    with preconditioning.

### Changed

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
