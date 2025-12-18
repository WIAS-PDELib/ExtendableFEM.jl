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
on a moving domain ``\Omega(t)`` (a sliding droplet) with Navier-slip and capillary boundary conditions,
see the reference below for details.

The ALE method computes a postprocessed displacement
``\mathbf{w}`` from ``\mathbf{u}`` by Laplace smoothing to prevent mesh folding. The smoothing
preserves the normal fluxes via a Lagrange multiplier.

The weak formulation characterizes ``(\mathbf{u},p,\mathbf{w},\lambda) \in V \times Q \times W \times \Lambda`` by
```math
\begin{aligned}
	(2\mu \varepsilon(\mathbf{u}), \varepsilon(\mathbf{v})) - (\mathrm{div} \mathbf{v}, p)\\
	+ \beta (\mathbf{u}, \mathbf{v})_{\Gamma_{solid}}
	+ \delta (\mathbf{u}, \mathbf{v})_{\Gamma_{contact}}
	& = \mathrm{Bo} (\mathbf{e}_x, \mathbf{v}) - (\nabla_s \mathbf{x}, \nabla_s \mathbf{v})_{\Gamma_{air}}
	+ \cos(\theta) (\nabla_s \mathbf{x}, \nabla_s \mathbf{v})_{\Gamma_{solid}}
	&& \quad \text{for all } \mathbf{v} \in V,\\
(\mathrm{div} \mathbf{u}, q) & = 0 && \quad \text{for all } q \in Q,\\
(\nabla \mathbf{w}, \nabla \mathbf{z}) + (\mathbf{z} \cdot \mathbf{n}, \lambda)_{\partial\Omega}
	& = 0
	&& \quad \text{for all } \mathbf{z} \in W,\\
    (\mathbf{w} \cdot \mathbf{n}, \psi)_{\partial\Omega}
	& = (\mathbf{u} \cdot \mathbf{n}, \psi)_{\partial\Omega}
	&& \quad \text{for all } \mathbf{\psi} \in \Lambda
\end{aligned}
```
Here, ``Bo``, ``\beta`` and ``\delta`` are dimensionless parameters
and ``\theta`` is the equilibration contact angle, see the reference below for details.

When the droplet speed, i.e. the mean velocity in x-direction,
becomes constant in time, a travelling wave solution is reached.
For the default parameters the result looks like this:

![](example295.png)


!!! reference

	"Resolving the microscopic hydrodynamics at the moving contact line",\
	A. K. Giri, P. Malgaretti, D. Peschka, and M. Sega,\
	Phys. Rev. Fluids 7 (2022),\
	[>Journal-Link<](https://doi.org/10.1103/PhysRevFluids.7.L102001)
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
    n = 2^(nref + 3)
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
        order = 2,          ## polynomial FEM order
        β = 1.0,            ## friction at liquid-solid interface
        δ = 1.0,            ## friction at contact line
        Bo = 0.5,           ## Bond number
        θ = π / 2,            ## equilibrium contact angle
        nsteps = 800,      ## number of ALE steps
        τ = 0.005,          ## ALE stepsize
        nrefs = 4,          ## mesh refinement level
        stationarity_target = 1.0e-2,
        Plotter = nothing, kwargs...
    )

    @info "δ = $δ, β = $(β)"

    ## grid
    xgrid, triple_nodes = initial_grid(nrefs)

    ## Stokes problem description
    PD = ProblemDescription("Stokes problem")
    u = Unknown("u"; name = "velocity")
    p = Unknown("p"; name = "pressure")
    x = Unknown("x"; name = "x")
    w = Unknown("w"; name = "extension")
    q = Unknown("q"; name = "Lagrange multiplier for normal flux")
    v = Unknown("v"; name = "comoving velocity")

    assign_unknown!(PD, u)
    assign_unknown!(PD, p)

    ## BLFs a(u,v), b(u,p)
    assign_operator!(PD, BilinearOperator([εV(u, 1.0)]; factor = 2, kwargs...))
    assign_operator!(PD, BilinearOperator([id(p)], [div(u)]; transposed_copy = 1, factor = -1, kwargs...))
    assign_operator!(PD, BilinearOperator([id(u)]; entities = ON_BFACES, regions = [2], factor = β, kwargs...))
    if abs(δ) > 0
        assign_operator!(PD, CallbackOperator(triple_junction_kernel!(triple_nodes, δ), [u]; kwargs...))
    end

    ## RHS
    assign_operator!(PD, LinearOperator(g!, [id(u)]; factor = Bo, kwargs...))
    assign_operator!(PD, LinearOperatorDG(surface_tension!, [grad(u)], [grad(x)]; entities = ON_BFACES, quadorder = 2, factor = -1, regions = [1], kwargs...))
    if abs(cos(θ)) > 1.0e-12
        assign_operator!(PD, LinearOperatorDG(surface_tension!, [grad(u)], [grad(x)]; entities = ON_BFACES, quadorder = 2, factor = cos(θ), regions = [2], kwargs...))
    end

    ## boundary conditions
    assign_operator!(PD, HomogeneousBoundaryData(u; regions = [2], mask = [0, 1, 1]))

    ## ALE problem description
    PDALE = ProblemDescription("ALE problem")
    assign_unknown!(PDALE, w)
    assign_unknown!(PDALE, q)
    assign_operator!(PDALE, BilinearOperator([grad(w)]; kwargs...))
    assign_operator!(PDALE, BilinearOperator([normalflux(w)], [id(q)]; transposed_copy = 1, regions = [1, 2], entities = ON_BFACES))
    assign_operator!(PDALE, LinearOperator([id(q)], [normalflux(u)]; regions = [1, 2], entities = ON_BFACES))
    assign_operator!(PDALE, HomogeneousBoundaryData(w; regions = [2], mask = [0, 1, 1]))

    ## prepare FESpace and solution vector
    FES = [FESpace{H1Pk{2, 2, order}}(xgrid), FESpace{H1Pk{1, 2, order - 1}}(xgrid), FESpace{H1Pk{1, 1, order - 1}, ON_BFACES}(xgrid), FESpace{H1P1{2}}(xgrid)]
    sol = FEVector([FES[1], FES[2]]; tags = [u, p])
    append!(sol, FES[1]; tag = w)
    append!(sol, FES[3]; tag = q)
    append!(sol, FES[4]; tag = x)

    SC = SolverConfiguration(PD, FES[[1, 2]]; init = sol, maxiterations = 1, verbosity = -1, timeroutputs = :none, kwargs...)
    SCALE = SolverConfiguration(PDALE, FES[[1, 3]]; init = sol, maxiterations = 1, verbosity = -1, timeroutputs = :none, kwargs...)

    ## prepare mean velocity integration
    vintegrate = ItemIntegrator([id(u, 1)]; piecewise = false)

    ## time loop
    time = 0.0
    mass = sum(xgrid[CellVolumes])
    v0, vprev = 0.0, 0.0
    for step in 1:nsteps
        time += τ

        ## redefine x for computation of the tangential identity grad(x)
        view(sol[x]) .= view(xgrid[Coordinates]', :)

        ## solve Stokes problem
        solve(PD, FES[[1, 2]], SC; kwargs...)

        ## solve ALE problem
        solve(PDALE, FES[[1, 3]], SCALE; kwargs...)

        ## displace mesh
        displace_mesh!(xgrid, sol[w]; magnify = τ)

        ## calculate droplet speed
        vprev = v0
        v0 = evaluate(vintegrate, sol)[1] / mass

        @info "STEP $step -------------
        time = $(Float16(time))
        droplet_speed = $(Float32(v0))
        stationarity = $(step > 1 ? Float32((vprev - v0) / τ) : String("init"))"

        if step > 1 && (vprev - v0 < τ * stationarity_target)
            @info "detected stationarity"
            break
        elseif step == nsteps
            @info "maximum number of ALE steps reached"
        end
    end

    ## compute comoving velocity (= u - v0)
    append!(sol, FES[1]; tag = v)
    view(sol[v]) .= view(sol[u])
    sol[v][1:FES[1].coffset] .-= v0

    plt = plot([grid(u), id(u), id(p), streamlines(v)], sol; Plotter = Plotter, levels = 0, colorbarticks = 7, rasterpoints = 21)

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
