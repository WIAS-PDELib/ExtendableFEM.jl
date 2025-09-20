#=

# 225 : Obstacle Problem Proximal Galerkin
([source code](@__SOURCE_URL__))

This example computes the solution ``u`` of the nonlinear obstacle problem that seeks the minimiser of the energy functional
```math
\begin{aligned}
	E(u) = \frac{1}{2} \int_\Omega \lvert \nabla u \rvert^2 dx - \int_\Omega f u dx
\end{aligned}
```
with some right-hand side ``f`` within the set of admissible functions that lie above an obstacle ``\chi``
```math
\begin{aligned}
	\mathcal{K} := \lbrace u \in H^1_0(\Omega) : u \geq \chi \rbrace.
\end{aligned}
```

Opposite to Example225 the solution is computed by the Proximal Galerkin method
that reformulates the problem into a solve of a series of nonlinear mixed problems.


!!! reference

	''Proximal Galerkin: A Structure-Preserving Finite Element Method for Pointwise Bound Constraints''
	Brendan Keith, Thomas M. Surowiec, Found Comput Math (2024)
	[>Link<](https://doi.org/10.1007/s10208-024-09681-8)


![](example227.png)
=#

module Example227_ObstacleProblemPG

using ExtendableFEM
using ExtendableGrids
using LinearAlgebra
using Test #hide

## define obstacle and penalty kernel
const b = 9 // 20
const d = sqrt(1 // 4 - b^2)
function χ(x)
    r = sqrt(x[1]^2 + x[2]^2)
    if r <= b
        return sqrt(1 // 4 - r^2)
    else
        return d + b^2 / d - b * r / d
    end
end

function R!(result, input, qpinfo)
    amphi1 = input - χ(qpinfo.x)
    return result[1] = (amphi1) * (ln(amphi1) - 1)
end
function ∇R!(result, input, qpinfo)
    return result[1] = χ(qpinfo.x) + exp(input[1])
end

function kernel_blf!(result, input, qpinfo)
    α = qpinfo.params[1]
    result .= α .* input
    return nothing
end

function kernel_lf!(result, input, qpinfo)
    result .= input
    return nothing
end

function main(; Plotter = nothing, nrefs = 6, α_0 = 1.0, order = 1, parallel = false, npart = 8, kwargs...)

    ## choose initial mesh
    xgrid = uniform_refine(grid_unitsquare(Triangle2D; scale = (2, 2), shift = (-0.5, -0.5)), nrefs)
    if parallel
        xgrid = partition(xgrid, RecursiveMetisPartitioning(npart = npart))
    end

    ## init proximal parameter
    params = [α_0]

    ## problem description
    PD = ProblemDescription()
    u = Unknown("u"; name = "u")
    ψ = Unknown("ψ"; name = "ψ")
    assign_unknown!(PD, u)
    assign_unknown!(PD, ψ)
    assign_operator!(PD, BilinearOperator(kernel_blf!, [grad(u)]; params = params, parallel = parallel, kwargs...))
    assign_operator!(PD, BilinearOperator([id(u)], [(id(ψ))]; transposed_copy = 1, store = true, parallel = parallel, kwargs...))
    assign_operator!(PD, NonlinearOperator(∇R!, [id(ψ)], [id(ψ)]; parallel = parallel, factor = -1, bonus_quadorder = 5, kwargs...))
    assign_operator!(PD, HomogeneousBoundaryData(u; regions = 1:4, kwargs...))

    ## create finite element space
    FES = [FESpace{H1Pk{1, 2, order}}(xgrid), FESpace{H1Pk{1, 2, order}}(xgrid)]

    ## add coupling matrix
    M = FEMatrix(FES[1], FES[2])
    b = FEVector(FES[1])
    assemble!(M, BilinearOperator([id(1)], [id(1)]))
    assign_operator!(PD, LinearOperator(b, [u]; factor = 1, kwargs...))

    ## solve
    sol = FEVector(FES; tags = PD.unknowns)
    sol_prev = FEVector(FES; tags = PD.unknowns)
    SC = nothing
    r = 3 // 2
    q = 3 // 2
    converged = false
    k = 0
    while !converged
        ## save previous solution and update right-hand side
        b.entries .= M.entries * view(sol[ψ])
        sol_prev.entries .= sol.entries

        ## solve nonlinear problem
        k += 1
        @info "Step $k: α = $(params[1])"
        sol, SC = solve(PD, FES, SC; init = sol, maxiterations = 20, return_config = true, kwargs...)
        dist = norm(view(sol[u]) .- view(sol_prev[u]))
        @info "dist = $dist"
        if dist < 1.0e-10
            converged = true
        else
            # increase proximal parameter
            params[1] = min(max(r^(q^k) - params[1]), 10^2)
        end
    end


    ## plot
    plt = plot([id(u), id(ψ)], sol; Plotter = Plotter, ncols = 3)

    return sol, plt
end

generateplots = ExtendableFEM.default_generateplots(Example227_ObstacleProblemPG, "example227.png") #hide
function runtests() #hide
    sol, plt = main(; μ = 1.0, nrefs = 2, order = 2) #hide
    @test maximum(sol.entries) ≈ 0.0033496680638875204 #hide
    return nothing #hide
end #hide
end # module
