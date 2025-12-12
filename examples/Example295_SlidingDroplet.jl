#=

# 295 : Sliding Droplet
([source code](@__SOURCE_URL__))

This problem describes the  motion of a liquid droplet sliding down a solid surface under gravity,
while maintaining its shape due to  surface tension and experiencing  slip at the solid boundary.
The problem is solved using an ALE (Arbitrary Lagrangian-Eulerian) approach where the mesh is moved
according to the extension velocity ``\mathbf{w}``.

It seeks a velocity ``\mathbf{u}``, a pressure ``p`` such that
```math
\begin{aligned}
	- 2\mu \varepsilon(\mathbf{u}) + \nabla p & = \mathbf{g} \\
	\mathrm{div} \mathbf{u} & = 0
\end{aligned}
```
on a moving domain ``\Omega(t)`` (a sliding droplet) with boundary conditions
```math
\begin{aligned}
	\mathbf{u} \cdot \mathbf{n} & = 0 && \quad \text{no penetration along liquid-solid interface} \\
	\mathbf{u} \cdot \mathbf{t} & = \mu_{sl}^{-1} \sigma_{nt} && \quad \text{Navier slip condition along liquid-solid interface} \\
	\sigma \mathbf{n} & = \gamma \kappa \mathbf{n} && \quad \text{surface tension at liquid-air interface} \\
	\sigma \mathbf{n} & = \gamma_{sl} \kappa \mathbf{n} && \quad \text{surface tension at liquid-solid interface}
\end{aligned}
```
where ``\sigma = -p I + 2\mu \varepsilon(\mathbf{u})`` is the stress tensor, ``\kappa`` is the curvature,
``\mathbf{n}`` is the normal vector and ``\mathbf{t}`` is the tangent vector.

Then, it computes the displacement ``\mathbf{w}`` from ``\mathbf{u}`` by Laplace smoothing to prevent mesh folding.

The weak formulation chracterizes ``(\mathbf{u},p,\mathbf{w}) \in V \times Q \times W`` by
```math
\begin{aligned}
	(2\mu \varepsilon(\mathbf{u}), \varepsilon(\mathbf{v})) + (\mathrm{div} \mathbf{v}, p)
	+ \gamma \tau (\nabla_s \mathbf{u}, \nabla_s \mathbf{v})_{\Gamma_{air}}
	+ \gamma_{sl} \tau (\nabla_s \mathbf{u}, \nabla_s \mathbf{v})_{\Gamma_{solid}}
	+ \mu_{sl} (\mathbf{u} \cdot \mathbf{t}, \mathbf{v} \cdot \mathbf{t})_{\Gamma_{solid}}
	& = (\mathbf{g}, \mathbf{v}) - \gamma (\nabla_s \mathbf{x}, \nabla_s \mathbf{v})_{\Gamma_{air}}
	- \gamma_{sl} (\nabla_s \mathbf{x}, \nabla_s \mathbf{v})_{\Gamma_{solid}}
	&& \quad \text{for all } \mathbf{v} \in V,\\
(\mathrm{div} \mathbf{u}, q) & = 0 && \quad \text{for all } q \in Q,\\
(\nabla \mathbf{w}, \nabla \mathbf{z}) + \epsilon^{-1} (\mathbf{w} \cdot \mathbf{n}, \mathbf{z} \cdot \mathbf{n})_{\partial\Omega}
	& = \epsilon^{-1} (\mathbf{u} \cdot \mathbf{n}, \mathbf{z} \cdot \mathbf{n})_{\partial\Omega}
	&& \quad \text{for all } \mathbf{z} \in W.
\end{aligned}
```

The computed solution for the default parameters looks like this:

![](example295.png)
=#

module Example295_SlidingDroplet

using ExtendableFEM
using ExtendableGrids
using ExtendableSparse
using LinearAlgebra
using Triangulate
using SimplexGridFactory
using GridVisualize
using Test #hide

## gavity
function g!(result, qpinfo)
    result[1] = 1.0
    return result[2] = 0.0
end

function surface_tension!(result, input, qpinfo)
    t1, t2 = qpinfo.normal[2], -qpinfo.normal[1] # tangent
    p1 = (input[1] * t1 + input[2] * t2)
    p2 = (input[3] * t1 + input[4] * t2)
    result[1] = p1 * t1
    result[2] = p1 * t2
    result[3] = p2 * t1
    return result[4] = p2 * t2
end

function initial_grid(nref; radius = 1)
    builder = SimplexGridBuilder(Generator = Triangulate)
    n = 2^(nref + 4)
    maxvol = 2.0^(-nref - 3)
    points = [point!(builder, radius * sin(t), radius * cos(t)) for t in range(-π / 2, π / 2, length = n)]

    facetregion!(builder, 1)
    for i in 1:(n - 1)
        facet!(builder, points[i], points[i + 1])
    end
    facetregion!(builder, 2)
    facet!(builder, points[end], points[1])

    return simplexgrid(builder, maxvolume = maxvol), [1, length(points)]
end

function main(;
        order = 2,
        g = 10,            ## gravity factor (in right direction)
        μ = 1,              ## bulk viscosity
        μ_sl = 0.1,         ## coefficient for Navier slip condition
        μ_dyn = 0,          ## dynamic contact angle
        γ_la = 1,              ## surface tension coefficient at liquid <> air interface
        γ_sl = 0,           ## surface tension coefficient at liquid <> solid interface
        nsteps = 70,       ## ALE steps
        τ = 1.0e-2,           ## ALE stepsize

        nrefs = 4,
        Plotter = nothing, kwargs...
    )

    ## grid
    xgrid, triple_nodes = initial_grid(nrefs)

    ## Stokes problem description
    PD = ProblemDescription("Stokes problem")
    u = Unknown("u"; name = "velocity")
    p = Unknown("p"; name = "pressure")
    x = Unknown("x"; name = "x")

    assign_unknown!(PD, u)
    assign_unknown!(PD, p)

    ## BLFs a(u,v), b(u,p)
    assign_operator!(PD, BilinearOperator([εV(u, 1.0)]; factor = 2 * μ, kwargs...))
    assign_operator!(PD, BilinearOperatorDG(surface_tension!, [grad(u)]; entities = ON_BFACES, factor = γ_la * τ, regions = [1], kwargs...))
    assign_operator!(PD, BilinearOperatorDG(surface_tension!, [grad(u)]; entities = ON_BFACES, factor = γ_sl * τ, regions = [2], kwargs...))
    assign_operator!(PD, BilinearOperator([id(p)], [div(u)]; transposed_copy = 1, factor = 1, kwargs...))
    assign_operator!(PD, BilinearOperator([id(u)]; entities = ON_BFACES, regions = [2], factor = μ_sl, kwargs...))
    if abs(μ_dyn) > 0
        assign_operator!(PD, CallbackOperator(ctriple_junction_kernel!(triple_nodes, μ_dyn), [u]; kwargs...))
    end

    ## RHS
    assign_operator!(PD, LinearOperator(g!, [id(u)]; factor = g, kwargs...))
    assign_operator!(PD, LinearOperatorDG(surface_tension!, [grad(u)], [grad(x)]; entities = ON_BFACES, quadorder = 2, factor = -γ_la, regions = [1], kwargs...))
    assign_operator!(PD, LinearOperatorDG(surface_tension!, [grad(u)], [grad(x)]; entities = ON_BFACES, quadorder = 2, factor = -γ_sl, regions = [2], kwargs...))

    ## boundary conditions
    assign_operator!(PD, HomogeneousBoundaryData(u; regions = [2], mask = [0, 1, 1]))

    ## ALE problem description
    penalty_normalfux = 1.0e5
    PDALE = ProblemDescription("ALE problem")
    w = Unknown("w"; name = "extension")
    assign_unknown!(PDALE, w)
    assign_operator!(PDALE, BilinearOperator([grad(w)]; kwargs...))
    assign_operator!(PDALE, BilinearOperator([normalflux(w)]; regions = [1, 2], entities = ON_BFACES, factor = penalty_normalfux))
    assign_operator!(PDALE, LinearOperator([normalflux(w)], [normalflux(u)]; regions = [1, 2], entities = ON_BFACES, factor = penalty_normalfux))
    assign_operator!(PDALE, HomogeneousBoundaryData(w; regions = [2], mask = [0, 1, 1]))

    ## prepare FESpace and solution vector
    FES = [FESpace{H1Pk{2, 2, order}}(xgrid), FESpace{H1Pk{1, 2, order - 1}}(xgrid)]
    sol = FEVector(FES; tags = [u, p])
    append!(sol, FES[1]; tag = w)
    append!(sol, FES[1]; tag = x)

    SC = SolverConfiguration(PD, FES[[1, 2]]; init = sol, maxiterations = 1, verbosity = -1, timeroutputs = :none, kwargs...)
    SCALE = SolverConfiguration(PDALE, FES[[1]]; init = sol, maxiterations = 1, verbosity = -1, timeroutputs = :none, kwargs...)

    ## prepare plot
    time = 0.0
    nodevals = nodevalues_view(sol[u])
    PE = PointEvaluator([id(u)], sol)
    plt = GridVisualizer(; Plotter = Plotter, layout = (3, 3), clear = true, resolution = (1600, 1200))
    gridplot!(plt[1, 1], xgrid; linewidth = 1, xlabel = "", ylabel = "", colorbar = false, title = "initial grid")
    dpt = (nsteps * τ) / 7
    times_to_plot = (1:7) * dpt

    ## time loop
    nextplot = 1
    for step in 1:nsteps
        @info "STEP $step, time = $(Float16(time))"

        ## compute x for computation of the tangential identity grad(x)
        interpolate!(sol[x], (result, qpinfo) -> (result .= qpinfo.x))

        ## solve Stokes problem
        solve(PD, FES[[1, 2]], SC; kwargs...)

        ## solve ALE problem
        solve(PDALE, FES[1], SCALE; kwargs...)

        ## displace mesh
        time += τ
        displace_mesh!(xgrid, sol[w]; magnify = τ)

        ## plot
        if time >= times_to_plot[nextplot]
            scalarplot!(plt[1 + Int(floor((nextplot) / 3)), 1 + (nextplot) % 3], xgrid, sqrt.(nodevals[1] .^ 2 .+ nodevals[2] .^ 2); levels = 7, title = "u (t = $(Float32(time)))")
            vectorplot!(plt[1 + Int(floor((nextplot) / 3)), 1 + (nextplot) % 3], xgrid, eval_func_bary(PE); vscale = 0.8, levels = 7, clear = false, title = "u (t = $(Float32(time)))")
            nextplot += 1
        end
    end

    ## calculate final droplet speed
    v0 = sum(sol.entries[1:FES[1].coffset]) / FES[1].coffset
    @info "droplet_speed = $v0"

    sol.entries[1:FES[1].coffset] .= sol.entries[1:FES[1].coffset] .- v0

    scalarplot!(plt[3, 3], xgrid, sqrt.(nodevals[1] .^ 2 .+ nodevals[2] .^ 2); levels = 7, title = "w (t = $(Float32(time)))")
    vectorplot!(plt[3, 3], xgrid, eval_func_bary(PE); vscale = 0.8, levels = 7, clear = false, title = "w (t = $(Float32(time)))")


    return sol, plt
end

## slip condition for triplet junctions (air-surface-solid)
function triple_junction_kernel!(triple_nodes, μ_dyn)
    factors = [1, 1, 0, 0]
    return function closure(A, b, args; assemble_matrix = true, kwargs...)
        FES = args[1].FES
        triple_dofs = copy(triple_nodes)
        append!(triple_dofs, triple_nodes .+ FES.coffset)
        if assemble_matrix
            for j in 1:length(triple_dofs)
                dof = triple_dofs[j]
                A[dof, dof] += factors[j] * μ_dyn
            end
            flush!(A)
        end
        return nothing
    end
end

generateplots = ExtendableFEM.default_generateplots(Example295_SlidingDroplet, "example295.png") #hide
end # module
