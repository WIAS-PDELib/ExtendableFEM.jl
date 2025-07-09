# Postprocessing and Visualization

A variety of postprocessing and visualization tools are available, both within this package and through related packages in the ecosystem. These tools enable you to analyze, visualize, and export finite element solutions and derived quantities.

## Related Packages

- [ExtendableGrids.jl](https://github.com/WIAS-PDELib/ExtendableGrids.jl):
    - Interface to [WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl) for exporting results to VTK format (for use with ParaView, VisIt, etc.).
    - `CellFinder` utility for locating the cell containing a given point.
- [ExtendableFEMBase.jl](https://github.com/WIAS-PDELib/ExtendableFEMBase.jl):
    - Node value interpolations.
    - `PointEvaluator` for evaluating solutions at arbitrary points in the domain.
    - `SegmentIntegrator` for integrating along 1D line intersections with the mesh.
    - Basic unicode plotting for quick inspection of results.
- [GridVisualize.jl](https://github.com/WIAS-PDELib/GridVisualize.jl):
    - Grid and scalar function plotting for simplicial grids in 1D, 2D, and 3D.
    - Supports multiple backends, e.g., PyPlot, GLMakie, PlutoVista.

## Plots and Tables

This package provides some convenient plotting and table-generation functions for visualizing and summarizing results:

```@autodocs
Modules = [ExtendableFEM]
Pages = ["plots.jl", "io.jl"]
Order   = [:type, :function]
```
