#=

# 260 : Axisymmetric Stokes
([source code](@__SOURCE_URL__))

This example solves the 3D stagnation point flow via 
the 2D axisymmetric formulation of the Navier--Stokes problem that seeks a velocity ``\mathbf{u} = (u_z, u_r)``
and pressure ``p`` such that
```math
\begin{aligned}
- \mu\left(\partial^2_r + r^{-1} \partial_r + \partial^2_z - r^{-2} \right) u_r
+ (u_r \partial_r + u_z \partial_z) u_r + \partial_r p & = \mathbf{f}_r\\
- \mu\left(\partial^2_r + r^{-1} \partial_r + \partial^2_z \right) u_z
+ (u_r \partial_r + u_z \partial_z) u_z + \partial_z p & = \mathbf{f}_z\\
(\partial_r + r^{-1})u_r + \partial_z u_z & = 0
\end{aligned}
```
with exterior force ``\mathbf{f}`` and some viscosity parameter ``\mu``.

The axisymmetric formulation assumes that the velocity in some
3D-domain, that is obtained by rotation of a 2D domain ``\Omega``,
only depends on the distance ``r`` to the rotation axis and the
``z``-coordinate tangential to the x-axis, but not on the angular coordinate
of the cylindric coordinates.
The implementation employs ``r``-dependent bilinear forms and a
Cartesian grid for the 2D ``(z,r)`` domain that is assumed to be rotated
around the ``r=0``-axis.

This leads to the weak formulation
```math
\begin{aligned}
a(u,v) + b(p,v) & = (f,v) \\
         b(q,u) & = 0
\end{aligned}
```
with the bilinear forms
```math
\begin{aligned}
a(u,v) := \int_{\Omega} \left( \nabla u : \nabla v + r^{-2} u_r v_r \right) r dr dz\\
b(q,v) := \int_{\Omega} q \left( \mathrm{div}(v) + r^{-1} u_r \right) r dr dz
\end{aligned}
```
where the usual Cartesian differential operators can be used. The factor ``2\pi`` from the
integral over the rotation angle drops out on both sides.

The computed solution for the default parameters looks like this:

![](example260.png)

=#

module Example260_AxisymmetricNavierStokesProblem

using ExtendableFEM
using ExtendableFEMBase
using ExtendableGrids
using ExtendableSparse
using SimplexGridFactory
using Triangulate
using Symbolics
using GridVisualize
using LaTeXStrings
using LinearAlgebra
using DataFrames
using CSV
using Test #hide


function H1BR_to_HDIVBDM1_axisymmetric_interpolator(FES_BR::FESpace{Tv, Ti, H1BR{2}}; axisymmetric_modifications = true, onlyRT0 = false) where {Tv, Ti}
    facedofs_V1::VariableTargetAdjacency{Int32} = FES_BR[FaceDofs]

    ## extract grid information
    xgrid = FES_BR.xgrid
    xFaceVolumes = xgrid[FaceVolumes]
    xCoordinates = xgrid[Coordinates]
    xFaceNodes = xgrid[FaceNodes]
    xFaceNormals = xgrid[FaceNormals]
    nfaces = num_sources(xFaceNormals)

    ## create reconstruction space
    FES_reconst = FESpace{HDIVBDM1{2}}(xgrid)

    ## create representation matrix of interpolation
    F = FEMatrix(FES_BR, FES_reconst)
    FE = F.entries
    w = zeros(Float64, 2)
    for face in 1:nfaces
        ## get radial coordinate
        w[1] = xCoordinates[1, xFaceNodes[1, face]]
        w[2] = xCoordinates[1, xFaceNodes[2, face]]
        wmid = (w[1] + w[2]) / 2

        ## RT0 dofs
        FE[facedofs_V1[1, face], face] = xFaceVolumes[face] * (1 // 6 * w[1] + 1 // 3 * wmid) * xFaceNormals[1, face]
        FE[facedofs_V1[2, face], face] = xFaceVolumes[face] * (1 // 6 * w[2] + 1 // 3 * wmid) * xFaceNormals[1, face]
        FE[facedofs_V1[3, face], face] = xFaceVolumes[face] * (1 // 6 * w[1] + 1 // 3 * wmid) * xFaceNormals[2, face]
        FE[facedofs_V1[4, face], face] = xFaceVolumes[face] * (1 // 6 * w[2] + 1 // 3 * wmid) * xFaceNormals[2, face]
        FE[facedofs_V1[5, face], face] = xFaceVolumes[face] * wmid

        ## BDM1 dofs
        if w[1] == 0 && axisymmetric_modifications
            for d in 1:5
                FE[facedofs_V1[d, face], nfaces + face] = 1 // 6 * FE[facedofs_V1[d, face], face]
            end
        elseif w[2] == 0 && axisymmetric_modifications
            for d in 1:5
                FE[facedofs_V1[d, face], nfaces + face] = -1 // 6 * FE[facedofs_V1[d, face], face]
            end
        elseif !onlyRT0
            FE[facedofs_V1[1, face], nfaces + face] = xFaceVolumes[face] * (-1 // 12 * w[1]) * xFaceNormals[1, face]
            FE[facedofs_V1[2, face], nfaces + face] = xFaceVolumes[face] * (1 // 12 * w[2]) * xFaceNormals[1, face]
            FE[facedofs_V1[3, face], nfaces + face] = xFaceVolumes[face] * (-1 // 12 * w[1]) * xFaceNormals[2, face]
            FE[facedofs_V1[4, face], nfaces + face] = xFaceVolumes[face] * (1 // 12 * w[2]) * xFaceNormals[2, face]
        end
    end
    flush!(FE)
    return F, FES_reconst
end


## exact data for problem by symbolics
function prepare_data(; μ = 1, vtype = 2, ptype = 1)

    @variables r z

    ## velocity u =(u_r, u_z)
    if vtype == 0
        u = [0 * r, 1]
    elseif vtype == 1
        u = [r, -2 * z]
    elseif vtype == 2
        u = [r^2, -3 * r * z]
    elseif vtype == 3
        u = [r * cos(pi * r / 2) * sin(pi * z + pi / 2), cos(pi * z + pi / 2) * ((2 / pi) * cos(pi * r / 2) - (r / 2) * sin(pi * r / 2))]
    end

    ## pressure
    if ptype == 0
        p = 0 * r
    elseif ptype == 1
        p = r^(7 / 4) - 2 / (15 / 4)
        p += z^2 - 1 / 3
    elseif ptype == 2
        p = log2(r)
    elseif ptype == 3
        s = 0.25
        p = r^s - 2 / (2 + s)
    elseif ptype == 4
        s = -0.5
        p = 10 * (r^s - 2 / (2 + s))
    elseif ptype == 5
        s = 3
        p = r^s - 2 / (2 + s)
    end

    @info "Computing axisymmetric Stokes data for u = $u and p = $p"

    ## gradient of velocity
    ∇u = Symbolics.jacobian(u, [r, z])
    ∇u_reshaped = [∇u[1, 1], ∇u[1, 2], ∇u[2, 1], ∇u[2, 2]]

    ## cylindrical Laplacian
    Δu = [
        (Symbolics.gradient(∇u[1, 1], [r]) + Symbolics.gradient(∇u[1, 2], [z]))[1] + ∇u[1, 1] / r,
        (Symbolics.gradient(∇u[2, 1], [r]) + Symbolics.gradient(∇u[2, 2], [z]))[1] + ∇u[2, 1] / r,
    ]

    ## right-hand side
    f = -μ * (Δu - [u[1] / r^2, 0]) + Symbolics.gradient(p, [r, z])
    @show f

    ## build functions
    p_eval = build_function(p, r, z, expression = Val{false})
    u_eval = build_function(u, r, z, expression = Val{false})
    ∇u_eval = build_function(∇u_reshaped, r, z, expression = Val{false})
    f_eval = build_function(f, r, z, expression = Val{false})

    return f_eval[2], u_eval[2], ∇u_eval[2], p_eval
end


function kernel_convection!(result, input, qpinfo)
    u, ∇u = view(input, 1:2), view(input, 3:6)
    r = qpinfo.x[1]
    result[1] = r * (∇u[1] * u[1] + ∇u[2] * u[2])
    result[2] = r * (∇u[3] * u[1] + ∇u[4] * u[2])
    return nothing
end

function kernel_stokes_axisymmetric!(result, u_ops, qpinfo)
    u, ∇u, p = view(u_ops, 1:2), view(u_ops, 3:6), view(u_ops, 7)
    r = qpinfo.x[1]
    μ = qpinfo.params[1]
    ## add Laplacian
    result[1] = μ / r * u[1] - p[1]
    result[2] = 0
    result[3] = μ * r * ∇u[1] - r * p[1]
    result[4] = μ * r * ∇u[2]
    result[5] = μ * r * ∇u[3]
    result[6] = μ * r * ∇u[4] - r * p[1]
    result[7] = -(r * (∇u[1] + ∇u[4]) + u[1])
    return nothing
end


function main(;
        domain = [(0, 1), (0, 1)],
        order = 1,
        vtype = 3,
        ptype = 3,
        μ = 0.1,
        levels = 1:6,
        nonlinear = false,
        uniform = false,
        quadorder = 10,
        quadorder_rhs = quadorder,
        quadorder_bnd = quadorder,
        quadorder_a = 0,
        operators = [:identity, :rt, :rt_plus, :bdm, :bdm_plus],
        Plotter = nothing, kwargs...
    )

    if isa(levels, AbstractArray)
        μ = μ * ones(length(levels))
        plot_vs = :level
    elseif isa(μ, AbstractArray)
        levels = levels * ones(Int, length(μ))
        plot_vs = :μ
    end

    @info μ, length(μ), levels

    ## unknowns
    u = Unknown("u"; name = "velocity")
    p = Unknown("p"; name = "pressure")
    ru = Unknown("ru"; name = "flux")
    nlevels = length(levels)
    NDofs = zeros(Int, nlevels)
    noperators = length(operators)
    Results = zeros(Float64, nlevels, 35)

    ## loop over reconstruction operators
    sol = nothing

    ## refinement loop
    for (l, level) in enumerate(levels)

        ## load problem data
        f_eval, u_eval, ∇u_eval, p_eval = prepare_data(; μ = μ[l], vtype = vtype, ptype = ptype)
        f!(result, qpinfo) = (f_eval(result, qpinfo.x[1], qpinfo.x[2]))
        p!(result, qpinfo) = (result[1] = p_eval(qpinfo.x[1], qpinfo.x[2]))
        u!(result, qpinfo) = (u_eval(result, qpinfo.x[1], qpinfo.x[2]))
        ∇u!(result, qpinfo) = (∇u_eval(result, qpinfo.x[1], qpinfo.x[2]))

        ## kernel for axisymmetric divergence
        function kernel_l2div(result, u_ops, qpinfo)
            idu, divu = view(u_ops, 1:2), view(u_ops, 3)
            result[1] = (qpinfo.x[1] * divu[1] + idu[1])^2 / qpinfo.x[1]
            return nothing
        end
        AxiDivIntegrator = ItemIntegrator(kernel_l2div, [id(u), div(u)]; quadorder = 10, resultdim = 1)

        ## kernel for pressure integration
        function kernel_pmean(result, ph, qpinfo)
            result[1] = qpinfo.x[1] * ph[1]
            return nothing
        end
        PMeanIntegrator = ItemIntegrator(kernel_pmean, [id(p)]; quadorder = 10, resultdim = 1)

        ## kernel for L2_{-1} norm of div(ru)
        function kernel_l2div_rv(result, div_ru, qpinfo)
            result[1] = div_ru[1]^2 / qpinfo.x[1]
            return nothing
        end
        DivIntegrator = ItemIntegrator(kernel_l2div_rv, [div(ru)]; quadorder = 10, resultdim = 1)

        ## kernel for L2_1 error of u
        function kernel_l2error(result, u_ops, qpinfo)
            u!(result, qpinfo)
            return result .= (result .- u_ops) .^ 2 * qpinfo.x[1]
        end
        L21ErrorIntegrator = ItemIntegrator(kernel_l2error, [id(u)]; quadorder = 10, kwargs...)

        ## kernel for L2_{-1} error of ru
        function kernel_l2m1error(result, u_ops, qpinfo)
            u!(result, qpinfo)
            return result .= (qpinfo.x[1] .* result .- u_ops) .^ 2 / qpinfo.x[1]
        end
        L2m1ErrorIntegrator = ItemIntegrator(kernel_l2m1error, [id(ru)]; quadorder = 10, kwargs...)

        ## kernel for L2 norm of ru on rotation axis
        function kernel_bnd_normalerror(result, u_ops, qpinfo)
            return result .= u_ops .^ 2
        end
        L2bndIntegrator = ItemIntegratorDG(kernel_bnd_normalerror, [this(id(ru))]; entities = ON_BFACES, regions = [4], resultdim = 2, quadorder = 10)

        ## kernel for L2_1 error of p
        function kernel_pl2error(result, ph, qpinfo)
            p!(result, qpinfo)
            return result .= (result .- ph) .^ 2 * qpinfo.x[1]
        end
        PErrorIntegratorExact = ItemIntegrator(kernel_pl2error, [id(p)]; quadorder = 10, kwargs...)

        ## kernel for energy error of u
        function kernel_energyerror(result, u_ops, qpinfo)
            u!(result, qpinfo)
            ∇u!(view(result, 3:6), qpinfo)
            result .-= u_ops
            result .= result .^ 2
            result[1] /= qpinfo.x[1]
            result[2] = 0
            result[3:6] .*= qpinfo.x[1]
            return nothing
        end
        EnergyErrorIntegratorExact = ItemIntegrator(kernel_energyerror, [id(u), grad(u)]; quadorder = 10, kwargs...)

        for (recid, reconstruction) in enumerate(operators)
            offset = (recid - 1) * 7

            ## grid
            if uniform
                hx = 2.0^-level
                hy = 2.0^-level
                xgrid = simplexgrid(domain[1][1]:hx:domain[1][2], domain[2][1]:hy:domain[2][2])
            else
                xgrid = simplexgrid(
                    Triangulate;
                    points = [domain[1][1] domain[2][1] ; domain[1][2] domain[2][1] ; domain[1][2] domain[2][2] ; domain[1][1] domain[2][2]]',
                    bfaces = [1 2 ; 2 3 ; 3 4 ; 4 1 ]',
                    bfaceregions = [1, 2, 3, 4],
                    regionpoints = [(domain[1][2] - domain[1][1]) / 2 (domain[2][2] - domain[2][1]) / 2;]',
                    regionnumbers = [1],
                    regionvolumes = [(domain[1][2] - domain[1][1]) * (domain[2][2] - domain[2][1]) * 4.0^(-level - 1)]
                )
            end

            ## FESpace
            if order == 1
                FES = [FESpace{H1BR{2}}(xgrid), FESpace{L2P0{1}}(xgrid)]
            else
                @error "only order = 1 currently implemented"
            end

            ## reconstruction space and operator
            F, FES_reconst = H1BR_to_HDIVBDM1_axisymmetric_interpolator(FES[1]; axisymmetric_modifications = reconstruction in [:bdm_plus, :rt_plus], onlyRT0 = reconstruction in [:rt, :rt_plus])

            ## problem description
            PD = ProblemDescription()
            assign_unknown!(PD, u)
            assign_unknown!(PD, p)
            assign_operator!(PD, BilinearOperator(kernel_stokes_axisymmetric!, [id(u), grad(u), id(p)]; bonus_quadorder = quadorder_a, params = [μ[l]], kwargs...)) #; jacobian = kernel_jacobian!))
            if nonlinear
                assign_operator!(PD, NonlinearOperator(kernel_convection!, [id(u)], [id(u), grad(u)]; bonus_quadorder = quadorder_a, kwargs...)) #; jacobian = kernel_jacobian!))
            end
            if reconstruction in [:rt, :rt_plus, :bdm, :bdm_plus]
                bR = FEVector(FES_reconst)
                assemble!(bR, LinearOperator(f!, [id(1)]; bonus_quadorder = quadorder_rhs, kwargs...); time = time)
                c = F.entries.cscmatrix * bR.entries
                assign_operator!(PD, LinearOperator(c, [u]; kwargs...))
            elseif reconstruction == :identity
                function fr!(result, qpinfo)
                    f!(result, qpinfo)
                    result .*= qpinfo.x[1]
                    return nothing
                end
                assign_operator!(PD, LinearOperator(fr!, [id(u)]; quadorder = quadorder_rhs, kwargs...))
            end
            assign_operator!(PD, HomogeneousBoundaryData(u; regions = [4], mask = (1, 0, 1)))
            assign_operator!(PD, InterpolateBoundaryData(u, u!; regions = [1, 2, 3], bonus_quadorder = quadorder_bnd))
            assign_restriction!(PD, BoundaryDataRestriction(u; regions = [4], value = 0, mask = (1, 0, 1)))
            assign_operator!(PD, FixDofs(p; dofs = [1]))

            ## solve
            sol = ExtendableFEM.solve(PD, FES; kwargs...)
            NDofs[l] = length(sol.entries)

            ## postprocess: append reconstruction of ru
            if reconstruction in [:bdm, :bdm_plus, :rt, :rt_plus]
                append!(sol, FES_reconst; tag = ru)
                view(sol[ru]) .= F.entries.cscmatrix' * view(sol[u])
            else
                FES_P3 = FESpace{H1Pk{2, 2, 3}}(xgrid)
                append!(sol, FES_P3; tag = ru)
                lazy_interpolate!(sol[ru], sol, [id(u)]; postprocess = (result, input, qpinfo) -> (result .= input .* qpinfo.x[1]))
            end

            ## correct pressure mean
            pmean = sum(evaluate(PMeanIntegrator, sol)) / sum(xgrid[CellVolumes])
            view(sol[p]) .-= 2 * pmean

            ## compute L2_1 error of velocity
            error = evaluate(L21ErrorIntegrator, sol)
            L2error = sqrt(sum(view(error, 1, :)) + sum(view(error, 2, :)))
            Results[l, 1 + offset] = L2error
            @info "||u - u_h||_{0,1} = $L2error"

            ## compute L2 error of velocity on rotation axis
            error = evaluate(L2bndIntegrator, sol)
            L2error_bnd = sqrt(sum(view(error, 1, :)) + sum(view(error, 2, :)))
            Results[l, 7 + offset] = L2error_bnd
            @info "||Π(ru_h)||_{0,Γ_rot} = $L2error_bnd"

            ## compute L2_{-1} error of velocity
            error = evaluate(L2m1ErrorIntegrator, sol)
            L2m1error = sqrt(sum(view(error, 1, :)) + sum(view(error, 2, :)))
            Results[l, 2 + offset] = L2m1error
            @info "||ru - Π(ru_h)||_{0,-1} = $L2m1error"

            ## compute energy error
            error = evaluate(EnergyErrorIntegratorExact, sol)
            Eerror = sqrt(sum(view(error, 1, :)) + sum(view(error, 2, :)) + sum(view(error, 3, :)) + sum(view(error, 4, :)) + sum(view(error, 5, :)) + sum(view(error, 6, :)))
            Results[l, 3 + offset] = Eerror
            @info "||u - u_h||_V = $Eerror"

            ## compute L2error of pressure
            error = evaluate(PErrorIntegratorExact, sol)
            pL2error = sqrt(sum(error))
            Results[l, 4 + offset] = pL2error
            @info "||p - p_h||_{0,1} = $pL2error"

            ## compute divergence in cylindrical coordinates by volume integrals
            div_error = sqrt(sum(evaluate(AxiDivIntegrator, sol)))
            Results[l, 5 + offset] = div_error
            @info "||div(u_h)|| = $div_error"

            ## compute divergence of reconstructed reduced velocity
            div_error_ru = sqrt(sum(evaluate(DivIntegrator, sol)))
            Results[l, 6 + offset] = div_error_ru
            @info "||div(Π(ru_h))||_{0,-1} = $div_error_ru"

            ## print results
            if plot_vs == :level
                @info "Results for reconstruction = $reconstruction at level = $level"
                ylabels = ["|| u - u_h ||_{0,1}", "|| ru - \\Pi (ru_h) ||_{0,-1}", "|| u - u_h ||_V", "|| p - p_h ||_{0,1}", "|| div u_h ||_{0,1}", "||div(Π(ru_h))||_{0,-1}"]
                print_convergencehistory(NDofs, Results[:, (offset + 1):(offset + 6)]; X_to_h = X -> X .^ (-1 / 2), ylabels = ylabels, xlabel = "ndof")
            end

            ## plot last solution
            # scalarplot!(plt[1, 1], id(u), sol; abs = true)
        end
    end


    ## produce plots
    if Plotter !== nothing

        ylabels = ["BR (Galerkin)", "BR (RT0)", "BR (modified RT0)", "BR (BDM1)", "BR (modified BDM1)"]

        ## plot convergence histories
        plt = GridVisualizer(; Plotter = Plotter, layout = (1, 1), clear = true, size = (600, 600))

        plot_convergencehistory!(
            plt[1, 1],
            plot_vs == :level ? NDofs : μ,
            Results[:, [3 + 7 * (m - 1) for m in 1:noperators]];
            xlabel = plot_vs == :level ? "ndof" : L"\mu",
            add_h_powers = plot_vs == :level ? [order] : [],
            X_to_h = X -> 2 * X .^ (-1 / 2),
            legend = :best,
            title = L"|| \mathbf{u} - \mathbf{u}_h ||_V",
            markershape = :star5,
            markersize = 20,
            linewidth = 3,
            fontsize = 32,
            # limits = (1e-3, 1e3),
            ylabels = ylabels,
        )
        scene = GridVisualize.reveal(plt)
        GridVisualize.save("energyerror_br_$(plot_vs == :μ ? "μvar" : μ[1])_vtype=$(vtype)_ptype=$(ptype)_uniform=$(uniform).png", scene; Plotter = Plotter)

        plot_convergencehistory!(
            plt[1, 1],
            plot_vs == :level ? NDofs : μ,
            Results[:, [1 + 7 * (m - 1) for m in 1:noperators]];
            xlabel = plot_vs == :level ? "ndof" : L"\mu",
            add_h_powers = plot_vs == :level ? [order + 1] : [],
            X_to_h = X -> X .^ (-1 / 2),
            legend = :best,
            title = L"|| \mathbf{u} - \mathbf{u}_h ||_{L^2_1}",
            marker = :star5,
            markersize = 20,
            linewidth = 3,
            fontsize = 32,
            # limits = (1e-6, 1e1),
            ylabels = ylabels,
        )
        scene = GridVisualize.reveal(plt)
        GridVisualize.save("l21error_br_$(plot_vs == :μ ? "μvar" : μ[1])_vtype=$(vtype)_ptype=$(ptype)_uniform=$(uniform).png", scene; Plotter = Plotter)

        plot_convergencehistory!(
            plt[1, 1],
            plot_vs == :level ? NDofs : μ,
            Results[:, [2 + 7 * (m - 1) for m in 1:noperators]];
            xlabel = plot_vs == :level ? "ndof" : L"\mu",
            add_h_powers = plot_vs == :level ? [order + 1] : [],
            X_to_h = X -> X .^ (-1 / 2),
            legend = :best,
            title = L"|| r\mathbf{u} - \Pi(r \mathbf{u}_h) ||_{L^2_{-1}}",
            marker = :star5,
            markersize = 20,
            linewidth = 3,
            fontsize = 32,
            # limits = (1e-6, 1e1),
            ylabels = ylabels,
        )
        scene = GridVisualize.reveal(plt)
        GridVisualize.save("l2m1error_br_$(plot_vs == :μ ? "μvar" : μ[1])_vtype=$(vtype)_ptype=$(ptype)_uniform=$(uniform).png", scene; Plotter = Plotter)

        plot_convergencehistory!(
            plt[1, 1],
            plot_vs == :level ? NDofs : μ,
            Results[:, [4 + 7 * (m - 1) for m in 1:noperators]];
            xlabel = plot_vs == :level ? "ndof" : L"\mu",
            add_h_powers = plot_vs == :level ? [order] : [],
            X_to_h = X -> 0.5 * X .^ (-1 / 2),
            legend = :best,
            title = L"|| p - p_h ||_{L^2_1}",
            marker = :star5,
            markersize = 20,
            linewidth = 3,
            fontsize = 32,
            # limits = (1e-4, 1e0),
            ylabels = ylabels,
        )
        scene = GridVisualize.reveal(plt)
        GridVisualize.save("perror_br_$(plot_vs == :μ ? "μvar" : μ[1])_vtype=$(vtype)_ptype=$(ptype)_uniform=$(uniform).png", scene; Plotter = Plotter)

        plot_convergencehistory!(
            plt[1, 1],
            plot_vs == :level ? NDofs : μ,
            Results[:, [7 + 7 * (m - 1) for m in 1:noperators]];
            xlabel = plot_vs == :level ? "ndof" : L"\mu",
            add_h_powers = plot_vs == :level ? [order] : [],
            X_to_h = X -> 0.5 * X .^ (-1 / 2),
            legend = :rc,
            title = L"|| Π(ru_h) ||_{L^2(\Gamma_\text{rot})}",
            marker = :star5,
            markersize = 20,
            linewidth = 3,
            fontsize = 32,
            # limits = (1e-4, 1e0),
            ylabels = ylabels,
        )
        scene = GridVisualize.reveal(plt)
        GridVisualize.save("rotaxiserror_br_$(plot_vs == :μ ? "μvar" : μ[1])_vtype=$(vtype)_ptype=$(ptype)_uniform=$(uniform).png", scene; Plotter = Plotter)

        plot([id(u), id(p), grid(u)], sol; Plotter = Plotter)
    end

    qlabels = [plot_vs == :μ ? "nu" : "ndofs", "L21 error", "L2m1 error", "H11 error", "L21 p error", "L2 div", "L2m1 div", "L2 bnd"]

    for (recid, reconstruction) in enumerate(operators)
        res = Results[:, 7 * (recid - 1) .+ (1:7)]
        if plot_vs == :level
            res = [NDofs res]
        else
            res = [μ res]
        end
        df = DataFrame(res, :auto)

        @info df
        CSV.write("results_$(reconstruction)_$(plot_vs == :μ ? "μvar" : μ[1])_vtype=$(vtype)_ptype=$(ptype)_uniform=$(uniform).csv", df, delim = ",", header = qlabels)
    end

    return Results, plt
end

generateplots = ExtendableFEM.default_generateplots(Example260_AxisymmetricNavierStokesProblem, "example260.png") #hide
function runtests() #hide
    errors, plt = main(; levels = 1:1) #hide
    @test all(Results .<= 1.0e-12) #hide
    return nothing #hide
end #hide
end # module
