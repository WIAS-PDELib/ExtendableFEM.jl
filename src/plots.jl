function GridVisualize.scalarplot!(
        p,
        op::Union{Tuple{Unknown, DataType}, Tuple{Int, DataType}},
        sol;
        abs = false,
        component = 1,
        title = typeof(op) <: Tuple{Unknown, DataType} ? String(op[1].identifier) : sol[op[1]].name,
        average_broken_plots = false,
        kwargs...
    )
    if !average_broken_plots && ExtendableFEMBase.broken(sol[op[1]].FES)
        broken_scalarplot!(p, sol[op[1]], op[2]; title, average_broken_plots, kwargs...)
    else
        return GridVisualize.scalarplot!(p, sol[op[1]].FES.dofgrid, view(nodevalues(sol[op[1]], op[2]; abs = abs), component, :); title = title, kwargs...)
    end
end

function GridVisualize.vectorplot!(p, op::Tuple{Unknown, DataType}, sol; title = String(op[1].identifier), kwargs...)
    return GridVisualize.vectorplot!(p, sol[op[1]].FES.dofgrid, eval_func_bary(PointEvaluator([op], sol)); title = title, kwargs...)
end
function GridVisualize.vectorplot!(p, op::Tuple{Int, DataType}, sol; title = sol[op[1]].name, kwargs...)
    return GridVisualize.vectorplot!(p, sol[op[1]].FES.dofgrid, eval_func_bary(PointEvaluator([op], sol)); title = title, kwargs...)
end

"""
````
function plot!(p::GridVisualizer, ops, sol; kwargs...)
````

Plots the operator evaluations ops of blocks in sol into the GridVisualizer.

"""
function plot!(
        p::GridVisualizer,
        ops,
        sol;
        keep = [],
        ncols = size(p.subplots, 2),
        reveal = true,
        kwargs...
    )
    col, row, id = 0, 1, 0
    for op in ops
        col += 1
        id += 1
        if col == ncols + 1
            col, row = 1, row + 1
        end
        while id in keep
            col += 1
            id += 1
            if col == ncols + 1
                col, row = 1, row + 1
            end
        end
        if typeof(op[1]) <: Unknown
            uname = op[2] == Identity ? String(op[1].identifier) : "$(op[2])(" * String(op[1].identifier) * ")"
        else
            uname = op[2] == Identity ? "$(sol[op[1]].name)" : "$(op[2])($(sol[op[1]].name))"
        end
        plot!(p, (row, col), op, sol; reveal = false, kwargs...)
    end
    if reveal
        GridVisualize.reveal(p)
    end
    return p
end


function plot!(
        p::GridVisualizer,
        index,
        op,
        sol;
        rasterpoints = 10,
        linewidth = 1,
        do_abs = true,
        do_vector_plots = true,
        title_add = "",
        reveal = true,
        average_broken_plots = false,
        kwargs...
    )
    if typeof(op[1]) <: Unknown
        uname = op[2] == Identity ? String(op[1].identifier) : "$(op[2])(" * String(op[1].identifier) * ")"
    else
        uname = op[2] == Identity ? "$(sol[op[1]].name)" : "$(op[2])($(sol[op[1]].name))"
    end
    if op[2] == "grid"
        gridplot!(p[index], sol[op[1]].FES.xgrid; linewidth = linewidth, title = uname * title_add, kwargs...)
    elseif op[2] == "dofgrid"
        gridplot!(p[index], sol[op[1]].FES.dofgrid; linewidth = linewidth, title = uname * title_add, kwargs...)
    elseif op[2] == "streamlines"
        PE = PointEvaluator([apply(op[1], Identity)], sol)
        streamplot!(p[index], sol[op[1]].FES.dofgrid, eval_func_bary(PE); rasterpoints = rasterpoints, title = uname * title_add, kwargs...)
    else
        ncomponents = get_ncomponents(sol[op[1]])
        edim = size(sol[op[1]].FES.xgrid[Coordinates], 1)
        resultdim = Length4Operator(op[2], edim, ncomponents)
        if resultdim == 1
            if !average_broken_plots && ExtendableFEMBase.broken(sol[op[1]].FES)
                ExtendableFEMBase.broken_scalarplot!(p[index], sol[op[1]], op[2]; title = uname * title_add, kwargs...)
            else
                GridVisualize.scalarplot!(p[index], sol[op[1]].FES.dofgrid, view(nodevalues(sol[op[1]], op[2]; abs = false), 1, :), title = uname * title_add; kwargs...)
            end
        elseif do_abs == true
            GridVisualize.scalarplot!(p[index], sol[op[1]].FES.dofgrid, view(nodevalues(sol[op[1]], op[2]; abs = true), 1, :), title = "|" * uname * "|" * title_add; kwargs...)
        else
            nv = nodevalues(sol[op[1]], op[2]; abs = false)
            for k in 1:resultdim
                if k > 1
                    col += 1
                    if col == ncols + 1
                        col, row = 1, row + 1
                    end
                end
                GridVisualize.scalarplot!(p[index], sol[op[1]].FES.dofgrid, view(nv, k, :), title = uname * " (component $k)" * title_add, kwargs...)
            end
        end
        if resultdim > 1 && do_vector_plots && do_abs == true && edim > 1
            GridVisualize.vectorplot!(p[index], sol[op[1]].FES.dofgrid, eval_func_bary(PointEvaluator([op], sol)); title = "|" * uname * "|" * " + quiver" * title_add, clear = false, kwargs...)
        end
    end
    if reveal
        GridVisualize.reveal(p)
    end
    return p
end

"""
````
function plot!(p::GridVisualizer, ops, sol; Plotter = nothing, kwargs...)
````

Plots the operator evaluations ops of blocks in sol with the specified Plotter module that is supported by GridVisualize (e.g. GLMakie, PyPlot, Plots)

"""
function plot(ops, sol; add = 0, Plotter = nothing, ncols = min(2, length(ops) + add), do_abs = true, width = (length(ops) + add) == 1 ? 400 : 800, height = 0, kwargs...)
    nplots = length(ops) + add
    for op in ops
        ncomponents = get_ncomponents(sol[op[1]])
        edim = size(sol[op[1]].FES.xgrid[Coordinates], 1)
        if !(op[2] in ["grid", "dofgrid", "streamlines"])
            resultdim = Length4Operator(op[2], edim, ncomponents)
            if resultdim > 1 && do_abs == false
                nplots += resultdim - 1
            end
        end
    end
    nrows = Int(ceil(nplots / ncols))
    if height == 0
        height = width / ncols * nrows
    end
    p = GridVisualizer(; Plotter = Plotter, layout = (nrows, ncols), clear = true, size = (width, height))
    return plot!(p, ops, sol; do_abs = do_abs, kwargs...)
end


function GridVisualize.vectorplot!(p, xgrid, op::Tuple{Union{Unknown, Int}, DataType}, sol; title = sol[op[1]].name, kwargs...)
    return GridVisualize.vectorplot!(p, xgrid, eval_func_bary(PointEvaluator([op], sol)); title = title, kwargs...)
end


"""
````
function plot_convergencehistory!(
	p::GridVisualizer, 
	X,
	Y;
	add_h_powers = [],
	X_to_h = X -> X,
	colors = [:blue, :green, :red, :magenta, :lightblue],
	title = "convergence history",
	legend = :best,
	ylabel = "",
	ylabels = [],
	xlabel = "ndofs",
	markershape = :circle,
	markevery = 1,
	clear = true,
	args...,
````

Plots a convergence history based on arrays X vs. Y into the GridVisualizer.

"""
function plot_convergencehistory!(
        target,
        X,
        Y;
        add_h_powers = [],
        X_to_h = X -> X,
        colors = [:blue, :green, :red, :magenta, :lightblue],
        title = "convergence history",
        legend = :best,
        ylabel = "",
        ylabels = [],
        xlabel = "ndofs",
        markershape = :circle,
        markevery = 1,
        clear = true,
        args...,
    )
    for j in 1:size(Y, 2)
        Xk = []
        Yk = []
        for k in 1:length(X)
            if Y[k, j] > 0
                push!(Xk, X[k])
                push!(Yk, Y[k, j])
            end
        end
        if length(ylabels) >= j
            label = ylabels[j]
        else
            label = "Data $j"
        end
        scalarplot!(
            target,
            simplexgrid(Xk),
            Yk;
            xlabel = xlabel,
            ylabel = ylabel,
            color = length(colors) >= j ? colors[j] : :black,
            clear = j == 1 ? clear : false,
            markershape = markershape,
            markevery = markevery,
            xscale = :log,
            yscale = :log,
            label = label,
            legend = legend,
            title = title,
            args...,
        )
    end
    for p in add_h_powers
        label = "h^$p"
        scalarplot!(target, simplexgrid(X), X_to_h(X) .^ p; linestyle = :dot, xlabel = xlabel, ylabel = ylabel, color = :black, clear = false, markershape = :none, xscale = :log, yscale = :log, label = label, legend = legend, title = title, args...)
    end
    return
end


"""
````
function plot_convergencehistory(X, Y; Plotter = nothing, kwargs...)
````

Plots a convergence history based on arrays X vs. Y into the GridVisualizer with the specified Plotter module (that needs to be supported by GridVisualize).

"""
function plot_convergencehistory(X, Y; Plotter = nothing, size = (800, 600), add_h_powers = [], X_to_h = X -> X, colors = [:blue, :green, :red, :magenta, :lightblue], legend = :best, ylabel = "", ylabels = [], xlabel = "ndofs", clear = true, args...)
    p = GridVisualizer(; Plotter = Plotter, layout = (1, 1), clear = true, size = size)
    return plot_convergencehistory!(p[1, 1], X, Y; add_h_powers = add_h_powers, X_to_h = X_to_h, colors = colors, legend = legend, ylabel = ylabel, ylabels = ylabels, xlabel = xlabel, clear = clear, args...)
end


function ExtendableFEMBase.nodevalues(op, sol; kwargs...)
    return nodevalues(sol[op[1]], op[2]; kwargs...)
end

## default function for generateplots for ExampleJuggler.jl
function default_generateplots(example_module, filename; kwargs...)
    return function closure(dir = pwd(); Plotter = nothing, kwargs...)
        ~, plt = example_module.main(; Plotter = Plotter, kwargs...)
        scene = GridVisualize.reveal(plt)
        return GridVisualize.save(joinpath(dir, filename), scene; Plotter = Plotter)
    end
end
