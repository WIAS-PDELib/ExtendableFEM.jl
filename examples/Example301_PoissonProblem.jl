#=

# 301 : Poisson-Problem
([source code](@__SOURCE_URL__))

This example computes the solution ``u`` of the two-dimensional Poisson problem
```math
\begin{aligned}
-\Delta u & = f \quad \text{in } \Omega
\end{aligned}
```
with right-hand side ``f(x,y) \equiv xy`` and homogeneous Dirichlet boundary conditions
on the unit cube domain ``\Omega`` on a given grid. The computed solution for the default
parameters looks like this:

![](example301.png)

This examples uses an iterative solver with an IncompleteLU preconditioner as the default solver.
It can be changed via the arguments 'method_linear' and 'precon_linear', see the runtests function
for more examples.

=#

module Example301_PoissonProblem

using ExtendableFEM
using ExtendableGrids
using LinearSolve
using IncompleteLU
using LinearAlgebra
using Test

function f!(result, qpinfo)
    result[1] = qpinfo.params[1]*(1.7^2 + 3.9^2)*sin(1.7*qpinfo.x[1])*cos(3.9*qpinfo.x[2])
    return nothing
end

function u!(result, qpinfo)
    result[1] = sin(1.7*qpinfo.x[1])*cos(3.9*qpinfo.x[2])
    return nothing
end

function main(;
    μ = 1.0,
    nrefs = 4, 
    method_linear = KrylovJL_GMRES(),
    precon_linear = method_linear == KrylovJL_GMRES() ? IncompleteLU.ilu : nothing,
    Plotter = nothing,
    kwargs...)

    ## problem description
    PD = ProblemDescription()
    u = Unknown("u"; name = "potential")
    assign_unknown!(PD, u)
    assign_operator!(PD, BilinearOperator([grad(u)]; factor = μ))
    assign_operator!(PD, LinearOperator(f!, [id(u)]; params = [μ]))
    assign_operator!(PD, InterpolateBoundaryData(u, u!; regions = 1:6))

    ## discretize
    xgrid = uniform_refine(grid_unitcube(Tetrahedron3D), nrefs)
    FES = FESpace{H1P2{1, 3}}(xgrid)

    ## solve
    sol = solve(PD, FES; method_linear, precon_linear, kwargs...)

    ## compute error
    function exact_error!(result, u, qpinfo)
        u!(result, qpinfo)
        result .-= u
        result .= result .^ 2
        return nothing
    end
    ErrorIntegratorExact = ItemIntegrator(exact_error!, [id(u)]; quadorder = 8)

    ## calculate error
    error = evaluate(ErrorIntegratorExact, sol)
    L2error = sqrt(sum(view(error, 1, :)))
    @info "L2 error = $L2error"

    ## plot
    plt = plot([id(u)], sol; Plotter = Plotter)

    return L2error, plt
end

generateplots = ExtendableFEM.default_generateplots(Example301_PoissonProblem, "example301.png")
function runtests()
    expected_error = 8.56e-5

    ## test direct solver
    L2error, plt = main(; nrefs = 4)
    @test L2error <= expected_error

    ## test iterative solver with IncompleteLU (fastest)
    method_linear = KrylovJL_GMRES()
    precon_linear = IncompleteLU.ilu
    L2error, plt = main(; method_linear, precon_linear, nrefs = 4)
    @test L2error <= expected_error

    ## test other iterative solvers
    method_linear = KrylovJL_GMRES(precs = (A, p) -> (Diagonal(A), I))
    precon_linear = nothing
    L2error, plt = main(; method_linear, precon_linear, nrefs = 4)
    @test L2error <= expected_error

    ## also working:
    ## method_linear = KrylovJL_GMRES(precs = (A, p) -> (AMGCLWrap.AMGPrecon(ruge_stuben(A)), I))
    ## method_linear = KrylovJL_GMRES(precs = (A, p) -> (AlgebraicMultigrid.aspreconditioner(ruge_stuben(A)), I))

    return nothing
end
end # module
