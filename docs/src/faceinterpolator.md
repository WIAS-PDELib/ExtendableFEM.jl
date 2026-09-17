# Face Interpolator

The face interpolator provides tools for evaluating and interpolating finite element functions on mesh faces. This is particularly useful for postprocessing tasks that require access to traces or jumps of functions across element boundaries.

## API Reference

```@autodocs
Modules = [ExtendableFEM]
Pages = ["common_operators/discface_interpolator.jl"]
Order   = [:type, :function]
```

## Example Usage (extracted from [Example210](module_examples/Example210_LshapeAdaptivePoissonProblem.md))

Suppose you want to compute the jumps of the gradient of a scalar-valued
Lagrange finite element function on the interior edges, e.g. to compute an a posteriori error estimator.

```julia
function gradnormalflux!(result, ∇u, qpinfo)
    result[1] = dot(∇u, qpinfo.normal)
end
NormalJumpProjector = FaceInterpolator(gradnormalflux!, [jump(grad(u))]; resultdim = 1, only_interior = true)
Jumps4Faces = evaluate!(NormalJumpProjector, sol)
```

See the [Example212](module_examples/Example212_PeriodicElasticity2D.md) for the complete example.
