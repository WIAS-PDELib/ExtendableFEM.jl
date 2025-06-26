#
# Print default dict for solver parameters into docstrings
#
function _myprint(dict::Dict{Symbol, Tuple{Any, String}})
    lines_out = IOBuffer()
    for k in sort!(collect(keys(dict)))
        v = dict[k]
        if typeof(v[1]) <: String
            println(lines_out, "  - `$(k)`: $(v[2]). Default: ''$(v[1])''\n")
        else
            println(lines_out, "  - `$(k)`: $(v[2]). Default: $(v[1])\n")
        end
    end
    return String(take!(lines_out))
end
#
# Update solver params from dict
#
function _update_params!(parameters, kwargs)
    for (k, v) in kwargs
        parameters[Symbol(k)] = v
    end
    return nothing
end

function center_string(S::String, L::Int = 8)
    if length(S) > L
        S = S[1:L]
    end
    while length(S) < L - 1
        S = " " * S * " "
    end
    if length(S) < L
        S = " " * S
    end
    return S
end


"""
$(TYPEDSIGNATURES)

Prints a formatted convergence history table for a sequence of data points, optionally in LaTeX format.

# Arguments
- `X`: Array of x-values (e.g., number of degrees of freedom, mesh sizes).
- `Y`: 2D array (or vector) of y-values (e.g., errors, residuals), with one row per entry in `X`.
- `X_to_h`: Function mapping `X` to a mesh size or similar for order calculation (default: identity).
- `ylabels`: Array of labels for each column of `Y` (default: auto-generated).
- `xlabel`: Label for the x-axis column (default: `"ndofs"`).
- `latex_mode`: If `true`, output is formatted as a LaTeX tabular environment (default: `false`).
- `separator`: Column separator for plain text or LaTeX (default: `"|"` or `"&"`).
- `order_seperator`: Separator between value and order columns (default: matches `separator`).
- `title`: (optional) Title for the table (default: none).
- `digits`: (optional) Number of significant digits for value formatting (default: `3`).

"""
function print_convergencehistory(
        X,
        Y;
        X_to_h = X -> X,
        ylabels = [],
        xlabel = "ndofs",
        latex_mode = false,
        separator = latex_mode ? "&" : "|",
        order_seperator = latex_mode ? "&" : "",
        title = "",
        digits = 3
    )
    # Input validation and error handling
    if ndims(Y) == 1
        Y = reshape(Y, :, 1)
    end
    if length(X) != size(Y, 1)
        error("Length of X ($(length(X))) must match number of rows in Y ($(size(Y, 1))).")
    end
    digits = max(1, Int(digits))
    # Automatic label generation
    if isempty(ylabels)
        ylabels = ["DATA $j" for j in 1:size(Y, 2)]
    end

    colwidth = max(maximum(length.(ylabels)), digits + 5)
    orderwidth = 6  # width for the order column
    xlabel = lpad(xlabel * " ", 12)
    # Print title if provided
    if !isempty(title)
        println(title)
    end
    if latex_mode
        tabular_argument = "c"
        for j in 1:size(Y, 2)
            tabular_argument *= "|cc"
        end
        @printf("\\begin{tabular}{%s}", tabular_argument)
    end
    @printf("\n%s%s", xlabel, separator)
    for j in 1:size(Y, 2)
        if length(ylabels) < j
            push!(ylabels, "DATA $j")
        end
        labelstr = lpad(ylabels[j], colwidth)
        orderstr = lpad("order", orderwidth)
        if j == size(Y, 2)
            @printf(" %s %s %s %s", labelstr, order_seperator, orderstr, latex_mode ? "" : separator)
        else
            @printf(" %s %s %s %s", labelstr, order_seperator, orderstr, separator)
        end
    end
    @printf("\n")
    if latex_mode
        @printf("\\\\\\hline")
    else
        @printf("%s%s", repeat("=", length(xlabel)), separator)
        for j in 1:size(Y, 2)
            @printf("%s%s", repeat("=", colwidth + orderwidth + 4), separator)
        end
    end
    @printf("\n")
    order = 0
    for j in 1:length(X)
        if eltype(X) <: Int
            @printf("   %7d  %s", X[j], separator)
        else
            @printf("  %.2e  %s", X[j], separator)
        end
        for k in 1:size(Y, 2)
            if j > 1
                order = -log(Y[j - 1, k] / Y[j, k]) / (log(X_to_h(X[j]) / X_to_h(X[j - 1])))
            end
            valstr = @sprintf("%*.*e", colwidth, digits, Y[j, k])
            if k == size(Y, 2)
                @printf(" %s%s %*.2f  %s", valstr, order_seperator, orderwidth, order, latex_mode ? "" : separator)
            else
                @printf(" %s%s %*.2f  %s", valstr, order_seperator, orderwidth, order, separator)
            end
        end
        if latex_mode
            @printf("\\\\")
        end
        @printf("\n")
    end
    if latex_mode
        @printf("\\end{tabular}")
    end
    return
end


"""
Prints a table with data X vs. Y

# Arguments
- `X`: Array of x-values (e.g., number of degrees of freedom, mesh sizes).
- `Y`: 2D array (or vector) of y-values (e.g., errors, residuals), with one row per entry in `X`.
- `ylabels`: Array of labels for each column of `Y` (default: auto-generated).
- `xlabel`: Label for the x-axis column (default: `"ndofs"`).
"""
function print_table(X, Y; ylabels = [], xlabel = "ndofs")
    xlabel = center_string(xlabel, 12)
    @printf("\n%s|", xlabel)
    for j in 1:size(Y, 2)
        if length(ylabels) < j
            push!(ylabels, "DATA $j")
        end
        @printf(" %s |", center_string(ylabels[j], 20))
    end
    @printf("\n")
    @printf("============|")
    for j in 1:size(Y, 2)
        @printf("======================|")
    end
    @printf("\n")
    for j in 1:length(X)
        if eltype(X) <: Int
            @printf("   %7d  |", X[j])
        else
            @printf("  %.2e  |", X[j])
        end
        for k in 1:size(Y, 2)
            @printf("    %.8e    |", Y[j, k])
        end
        @printf("\n")
    end
    return
end
