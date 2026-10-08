module AdaptivePressureRobustStokes

using Dates
using ExtendableFEM
using ExtendableGrids
using GridVisualize
using LaTeXStrings
using LinearAlgebra
using Printf
using PythonPlot
using Symbolics

default_plotter!(PythonPlot)

struct StokesData
    exact_u!::Function
    exact_∇u!::Function
    exact_p!::Function
    rhs!::Function
    exact_curl_f!::Function
end

function map_expression_to_eval_func(
        expr,
        variables,
        transformation::Union{Function, Nothing} = nothing
    )::Function
    compiled = build_function(expr, variables..., expression = Val{false})
    eval_func! = isa(compiled, Tuple) ? compiled[2] : (result, values...) -> (result[1] = compiled(values...))
    get_values = !isnothing(transformation) ? transformation : (qpinfo) -> qpinfo.x
    return (result, qpinfo) -> eval_func!(result, get_values(qpinfo)...)
end

function prepare_data_for_example_1(; ν::Real)::StokesData
    @variables x y

    ξ = x^2 * (x - 1)^2 * y^2 * (y - 1)^2
    ∇ξ = Symbolics.gradient(ξ, [x, y])
    u = [-∇ξ[2], ∇ξ[1]]
    p = x^5 + y^5 - 1 / 3
    ∇u = Symbolics.jacobian(u, [x, y])
    ∇u_reshaped = [∇u[1, 1], ∇u[1, 2], ∇u[2, 1], ∇u[2, 2]]
    Δu = [
        Symbolics.derivative(∇u[1, 1], x) + Symbolics.derivative(∇u[1, 2], y),
        Symbolics.derivative(∇u[2, 1], x) + Symbolics.derivative(∇u[2, 2], y),
    ]
    ∇p = Symbolics.gradient(p, [x, y])
    f = -ν * Δu + ∇p
    curl_f = Symbolics.derivative(f[2], x) - Symbolics.derivative(f[1], y)

    return StokesData(
        map(
            (expr) -> map_expression_to_eval_func(expr, [x, y]),
            [u, ∇u_reshaped, p, f, curl_f],
        )...,
    )
end

function prepare_data_for_example_2(; ν::Real)::StokesData
    @variables ϕ r

    α = 856399 / 1572864
    ω = 3 * pi / 2
    ψ = (1 / (α + 1)) * sin((α + 1) * ϕ) * cos(α * ω) - cos((α + 1) * ϕ)
    -(1 / (α - 1)) * sin((α - 1) * ϕ) * cos(α * ω) + cos((α - 1) * ϕ)
    dψ_dϕ = Symbolics.derivative(ψ, ϕ)
    d3ψ_dϕ3 = Symbolics.derivative(Symbolics.derivative(dψ_dϕ, ϕ), ϕ)
    p0 = -ν * r^(α - 1) * ((1 + α)^2 * dψ_dϕ + d3ψ_dϕ3) / (1 - α)
    u = [
        r^α * ((α + 1) * sin(ϕ) * ψ + cos(ϕ) * dψ_dϕ),
        r^α * (sin(ϕ) * dψ_dϕ - (α + 1) * cos(ϕ) * ψ),
    ]
    ∇u_rϕ = Symbolics.jacobian(u, [r, ϕ])
    ∇u = [
        ∇u_rϕ[1, 1] * cos(ϕ) - ∇u_rϕ[1, 2] * sin(ϕ) / r ∇u_rϕ[1, 1] * sin(ϕ) + ∇u_rϕ[1, 2] * cos(ϕ) / r;
        ∇u_rϕ[2, 1] * cos(ϕ) - ∇u_rϕ[2, 2] * sin(ϕ) / r ∇u_rϕ[2, 1] * sin(ϕ) + ∇u_rϕ[2, 2] * cos(ϕ) / r
    ]
    ∇u_reshaped = [∇u[1, 1], ∇u[1, 2], ∇u[2, 1], ∇u[2, 2]]
    p_plus = sin(r^2 * cos(ϕ) * sin(ϕ) * pi)
    ∇p_plus_rϕ = Symbolics.gradient(p_plus, [r, ϕ])
    f = [
        ∇p_plus_rϕ[1] * cos(ϕ) - ∇p_plus_rϕ[2] * sin(ϕ) / r,
        ∇p_plus_rϕ[1] * sin(ϕ) + ∇p_plus_rϕ[2] * cos(ϕ) / r,
    ]
    curl_f = 0 * r
    p = p0 + p_plus - 1 / 3 * 0.52466306757531903964

    function map_cartesian_to_polar(qpinfo)
        x, y = qpinfo.x
        return sqrt(x^2 + y^2), (2 * pi + atan(y, x)) % (2 * pi)
    end

    return StokesData(
        map(
            (expr) -> map_expression_to_eval_func(expr, [r, ϕ], map_cartesian_to_polar),
            [u, ∇u_reshaped, p, f, curl_f],
        )...,
    )
end

struct ExampleConfiguration
    prepare_data::Union{typeof(prepare_data_for_example_1), typeof(prepare_data_for_example_2)}
    xgrid::ExtendableGrid{Float64, Int32}
    boundary_regions::UnitRange{Int}
    options::Dict
end

example_1 = ExampleConfiguration(
    prepare_data_for_example_1,
    uniform_refine(grid_unitsquare(Triangle2D), 3),
    1:4,
    Dict("example_name" => "Example 1")
)

example_2 = ExampleConfiguration(
    prepare_data_for_example_2,
    uniform_refine(grid_lshape(Triangle2D), 3),
    1:8,
    Dict("example_name" => "Example 2")
)

# Kernel function for face interpolation of normal jumps of gradient
function gradnormalflux_kernel!(result, ∇u, qpinfo)
    result .= [∇u[1] ∇u[2]; ∇u[3] ∇u[4]] * qpinfo.normal
    return nothing
end

# Kernel function for face interpolation of tangental jumps of the Laplacian
function laplacetangentflux_kernel!(result, Δu, qpinfo)
    result .= Δu[1] * qpinfo.normal[2] - Δu[2] * qpinfo.normal[1]
    return nothing
end

function η_jump2_kernel!(result, laplacejump, qpinfo)
    result .= qpinfo.volume^3 * laplacejump .^ 2
    return nothing
end

function η_jump_kernel!(result, gradjump, qpinfo)
    result .= qpinfo.volume * sum(gradjump .^ 2)
    return nothing
end

function η_cons1_kernel!(result, Δu, qpinfo)
    result .= qpinfo.volume .* sum(Δu .^ 2)
    return nothing
end

function run_single(;
        configuration::ExampleConfiguration = example_1,
        ν::Real = 10^-3,
        reconstruct::Bool = true,
        pressure_robust_estimate::Bool = true,
        maxdofs::Int = 10000,
        maxlevel::Int = typemax(Int64),
        bonus_quadorder::Int = 5,
        Plotter = default_plotter(),
        kwargs...,
    )
    @assert ν > 0
    # Data preparation
    (; prepare_data, xgrid, boundary_regions, options) = configuration
    (; exact_u!, exact_∇u!, exact_p!, rhs!, exact_curl_f!) = prepare_data(; ν)

    # Kernel function for the Stokes problem
    function stokes_kernel!(result, input, _)
        result[1] = ν * input[1] - input[5]
        result[2] = ν * input[2]
        result[3] = ν * input[3]
        result[4] = ν * input[4] - input[5]
        result[5] = -(input[1] + input[4])
        return nothing
    end

    # Kernel functions for error calculations
    function velocity_error!(result, u, qpinfo)
        exact_u!(view(result, 1:2), qpinfo)
        exact_∇u!(view(result, 3:6), qpinfo)
        result .-= u
        result .= result .^ 2
        return nothing
    end

    function pressure_error!(result, p, qpinfo)
        exact_p!(view(result, 1), qpinfo)
        result .-= p
        result .= result .^ 2
        return nothing
    end

    # Kernel function for calculation of η_curl a posteriori error term
    function η_curl_kernel!(result, _, qpinfo)
        exact_curl_f!(view(result, 1), qpinfo)
        result .= qpinfo.volume^2 * result .^ 2
        return nothing
    end

    function η_vol_kernel!(result, input, qpinfo)
        Laplace_u = view(input, 1:2)
        ∇q = view(input, 3:4)
        f = Array{Float64}(undef, 2)
        rhs!(f, qpinfo)
        result .= qpinfo.volume * sum((f - ∇q + ν * Laplace_u) .^ 2)
        return nothing
    end

    # for classical problem (reconstruction = false) together with the pressure-robust a posteriori estimate
    function η_cons3_kernel!(result, input, qpinfo)
        Laplace_u = view(input, 1:2)
        f = Array{Float64}(undef, 2)
        rhs!(f, qpinfo)
        result .= qpinfo.volume * sum((f + ν * Laplace_u) .^ 2)
        return nothing
    end

    # define interpolators and item integrators for error estimation and calculation
    u = Unknown("u", name = "velocity", symbol_ansatz = "u", symbol_test = "v", dim = 2)
    p = Unknown("p", name = "pressure", symbol_ansatz = "p", symbol_test = "q", dim = 1)
    Δu = Unknown("Δu", name = "Laplacian of velocity", dim = 1)

    exact_velocity_error_integrator = ItemIntegrator(
        velocity_error!, [id(u), grad(u)];
        quadorder = bonus_quadorder, kwargs...,
    )
    exact_pressure_error_integrator = ItemIntegrator(
        pressure_error!, [id(p)];
        quadorder = bonus_quadorder, kwargs...,
    )
    η_jump_integrator = ItemIntegrator(
        η_jump_kernel!, [id(1)];
        resultdim = 1, quadorder = bonus_quadorder, entities = ON_FACES, kwargs...,
    )
    normal_jump_projector = FaceInterpolator(
        gradnormalflux_kernel!, [jump(grad(u))];
        resultdim = 2, quadorder = bonus_quadorder, only_interior = true, kwargs...,
    )
    tangent_jump_projector = FaceInterpolator(
        laplacetangentflux_kernel!, [jump(id(Δu))];
        resultdim = 1, quadorder = bonus_quadorder, entities = ON_FACES, kwargs...,
    )
    η_curl_integrator = ItemIntegrator(
        η_curl_kernel!, [Δ(u)];
        resultdim = 1, quadorder = bonus_quadorder, kwargs...,
    )
    η_jump2_integrator = ItemIntegrator(
        η_jump2_kernel!, [id(1)];
        resultdim = 1, quadorder = bonus_quadorder, entities = ON_FACES, kwargs...,
    )
    η_cons0_integrator = L2NormIntegrator(
        [div(u)];
        quadorder = bonus_quadorder, kwargs...,
    )
    η_cons1_integrator = ItemIntegrator(
        η_cons1_kernel!, [id(Δu)];
        resultdim = 1, quadorder = bonus_quadorder,
    )
    η_cons3_integrator = ItemIntegrator(
        η_cons3_kernel!, [id(Δu)];
        resultdim = 1, quadorder = bonus_quadorder,
    )
    η_vol_integrator = ItemIntegrator(
        η_vol_kernel!, [id(Δu), grad(p)];
        resultdim = 1, quadorder = bonus_quadorder,
    )

    FETypes = (H1BR{2}, L2P0{1}) # use BR finite element pair

    # initialize data structures for loop
    η_results = Array{Float64}(undef, 0)
    u_l2error_results = Array{Float64}(undef, 0)
    ∇u_l2error_results = Array{Float64}(undef, 0)
    p_l2error_results = Array{Float64}(undef, 0)
    ndofs = Array{Int64}(undef, 0)
    sol = nothing

    # SOLVE -> ESTIMATE -> MARK -> REFINE loop
    for level in Iterators.countfrom(1)
        # Create Stokes problem description
        # refined grid -> more boundary data -> recreate problem description
        PD = ProblemDescription("Stokes problem")
        assign_unknown!(PD, u)
        assign_unknown!(PD, p)

        id_u = reconstruct ? apply(u, Reconstruct{HDIVRT0{2}, Identity}) : id(u)
        assign_operator!(PD, BilinearOperator(stokes_kernel!, [grad(u), id(p)]))
        assign_operator!(PD, LinearOperator(rhs!, [id_u]; bonus_quadorder))
        assign_operator!(
            PD,
            InterpolateBoundaryData(u, exact_u!, ON_BFACES; regions = boundary_regions, bonus_quadorder)
        )
        assign_restriction!(PD, ZeroMeanValueRestriction(p))

        # Solve the Stokes problem
        @time begin
            FES = [FESpace{FETypes[1]}(xgrid), FESpace{FETypes[2]}(xgrid)]
            sol = solve(PD, FES)
            push!(ndofs, length(sol[1]) + length(sol[2]))
            println("\t ndof = ", ndofs[level])
            print("@time solver =")
        end

        # Calculate error statistics
        error = evaluate(exact_velocity_error_integrator, sol)
        push!(u_l2error_results, sqrt(sum(view(error, 1, :)) + sum(view(error, 2, :))))
        push!(∇u_l2error_results, sqrt(sum(sum(view(error, 3:6, :)))))
        evaluate!(error, exact_pressure_error_integrator, sol)
        push!(p_l2error_results, sqrt(sum(view(error, 1, :))))

        append!(sol, FESpace{L2P0{2}}(xgrid); tag = Δu)
        lazy_interpolate!(sol[Δu], sol, [Δ(u)]) # calculation of Δu ≡ (div ∘ ∇)u

        # Estimate the a posteriori error bound
        # η_curl = 0, as the Bernardi-Raugel element is used for the velocity space.
        # η_avg = 0, since σ = ∇u_h is chosen in the a posteriori estimator.
        # η_cons2 = 0, by choice of the Bernardi-Raugel reconstruction operator.
        normal_jumps_for_faces = evaluate!(normal_jump_projector, sol)
        η_jump = evaluate(η_jump_integrator, normal_jumps_for_faces)
        η_cons0 = evaluate(η_cons0_integrator, sol)
        η = η_jump
        η_cells = η_cons0
        if pressure_robust_estimate
            η_curl = evaluate(η_curl_integrator, sol)
            η_cells += η_curl
            if reconstruct
                tangent_jumps_for_faces = evaluate!(tangent_jump_projector, sol)
                η_jump2 = evaluate(η_jump2_integrator, tangent_jumps_for_faces)
                η += η_jump2
                η_cons1 = evaluate(η_cons1_integrator, sol)
                η_cells += η_cons1
            else
                η_cons3 = evaluate(η_cons3_integrator, sol)
                η_cells += η_cons3
            end
        else
            η_vol = evaluate(η_vol_integrator, sol)
            η_cells += ν^-2 * η_vol
            if reconstruct
                η_cons1 = evaluate(η_cons1_integrator, sol)
                η_cells += η_cons1
            end
        end

        for (index, face) in enumerate(eachcol(xgrid[FaceCells]))
            η[index] += η_cells[face[1]]
            if face[2] > 0
                η[index] += η_cells[face[2]]
            end
        end

        push!(η_results, sqrt(sum(η)))
        println(
            "\t    η =  $(η_results[level])" *
                "\n\t    ||u - u_h|| =  $(u_l2error_results[level])" *
                "\n\t    ||∇(u - u_h)|| =  $(∇u_l2error_results[level])" *
                "\n\t    ||p - p_h|| =  $(p_l2error_results[level])"
        )

        if ndofs[level] >= maxdofs || level >= maxlevel ||
                (level >= 2 && ndofs[level - 1] == ndofs[level])
            break
        end

        # Mark and refine steps
        @time begin
            facemarker = bulk_mark(xgrid, view(η, :), 0.5; indicator_AT = ON_FACES)
            xgrid = RGB_refine(xgrid, facemarker)
            print("@time refine =")
        end
    end

    # Create plots
    example_name = get(options, "example_name", nothing)
    estimator_type = pressure_robust_estimate ? "'pressure robust'" : "'classical'"
    title_details = !isnothing(example_name) ? (
            Printf.@sprintf(
                " for %s with ν = %.e, reconstruction=%s, a posteriori estimate=%s",
                example_name, ν, string(reconstruct), estimator_type
            )
        ) : ""
    solution_plt = ExtendableFEM.plot([id(u), id(p)], sol; Plotter)
    grid_plt = GridVisualizer(layout = (1, 1), size = (1500, 840))
    ExtendableFEM.plot!(
        grid_plt, [ExtendableFEM.grid(u)], sol;
        Plotter, title = "Grid$title_details, ndofs=$(ndofs[end])", show_colorbar = false
    )
    convergence_history_plt = GridVisualizer(layout = (1, 1), size = (1500, 840))
    ExtendableFEM.plot_convergencehistory!(
        convergence_history_plt,
        ndofs,
        [u_l2error_results ∇u_l2error_results p_l2error_results η_results];
        add_h_powers = [1, 2],
        X_to_h = X -> X .^ (-1 / 2),
        title = "Convergence history$title_details",
        xlabel = "number of degrees of freedom",
        ylabel = "value",
        ylabels = [
            "\$\\|u - u_h\\|_{L^2}\$",
            "\$\\|\\nabla(u - u_h)\\|_{L^2}\$",
            "\$\\|p - p_h\\|_{L^2}\$",
            "\$\\eta\$",
        ],
    )

    return sol, solution_plt, grid_plt, convergence_history_plt
end

function run_multiple(;
        configurations::Vector{ExampleConfiguration} = [example_1, example_2],
        Ν::Vector{Float64} = [10.0^-k for k in 0:6],
        reconstruct::Bool = true,
        pressure_robust_estimate::Bool = true,
        maxdofs::Int = 100000,
        maxlevel::Int = typemax(Int64),
        base_path = replace(string(Dates.now()), ":" => "_")
    )
    if !isdir(base_path)
        mkdir(base_path)
    end

    for configuration in configurations
        (; options) = configuration
        example_path = joinpath(base_path, replace(lowercase(options["example_name"]), " " => "_"))
        if !isdir(example_path)
            mkdir(example_path)
        end

        for ν in Ν
            single_run_path = joinpath(example_path, string(ν))
            if !isdir(single_run_path)
                mkdir(single_run_path)
            end

            _, solution_plt, grid_plt, convergence_history_plt =
                run_single(; configuration, ν, maxdofs, maxlevel, reconstruct, pressure_robust_estimate)

            GridVisualize.save(joinpath(single_run_path, "solution.svg"), solution_plt)
            GridVisualize.save(joinpath(single_run_path, "grid.svg"), grid_plt)
            GridVisualize.save(
                joinpath(single_run_path, "convergence_history.svg"),
                convergence_history_plt
            )
        end
    end

    return nothing
end

end
