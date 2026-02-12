default_statistics(Tv = Float64, Ti = Int64) = Dict{Symbol, Any}(
    :timeroutputs => TimerOutput(),
    :linear_residuals => Tv[],
    :nonlinear_residuals => Tv[],
    :matrix_nnz => Ti[],
    :restriction_residuals => Vector{Tv}[]
)

mutable struct SolverConfiguration{AT <: AbstractMatrix, bT, xT}
    PD::ProblemDescription
    A::AT ## stores system matrix
    b::bT ## stores right-hand side
    sol::xT ## stores solution
    tempsol::xT ## temporary solution
    res::xT
    freedofs::Vector{Int} ## stores indices of free dofs
    statistics::Dict{Symbol, Any}
    linsolver::Any
    unknown_ids_in_sol::Array{Int, 1}
    unknowns::Array{Unknown, 1}
    unknowns_reduced::Array{Unknown, 1}
    offsets::Array{Int, 1}  ## offset for each unknown that is solved
    parameters::Dict{Symbol, Any} # dictionary with user parameters
end

"""
````
residuals(S::SolverConfiguration)
````

returns the vector with the residuals of all iterations
"""
residuals(S::SolverConfiguration) = S.statistics[:nonlinear_residuals]

"""
````
residual(S::SolverConfiguration)
````

returns the residual of the last solve
"""
residual(S::SolverConfiguration) = S.statistics[:nonlinear_residuals][end]


"""
````
timeroutputs(S::SolverConfiguration)
````

returns TimerOutputs object that contains detailed information on solving and assembly times
"""
timeroutputs(S::SolverConfiguration) = S.statistics[:timeroutputs]


"""
````
lastmatrix(S::SolverConfiguration)
````

returns the currently stored system matrix
"""
lastmatrix(S::SolverConfiguration) = S.A

"""
````
lastrhs(S::SolverConfiguration)
````

returns the currently stored right-hand side
"""
lastrhs(S::SolverConfiguration) = S.b


#
# Default context information with help info.
#
default_solver_kwargs() = Dict{Symbol, Tuple{Any, String}}(
    :abstol => (1.0e-11, "abstol for linear solver (if iterative)"),
    :check_matrix => (false, "check matrix for symmetry and positive definiteness and largest/smallest eigenvalues"),
    :compress_restrictions => (true, "compress the column space of the given restriction matrices via QR to avoid overdetermined systems"),
    :constant_matrix => (false, "matrix is constant (skips reassembly and refactorization in solver)"),
    :constant_rhs => (false, "right-hand side is constant (skips reassembly)"),
    :damping => (0, "amount of damping, value should be between in (0,1)"),
    :inactive => (Array{Unknown, 1}([]), "inactive unknowns (are made available in assembly, but not updated in solve)"),
    :init => (nothing, "initial solution (also used to save the new solution)"),
    :initialized => (false, "linear system in solver configuration is already assembled (turns true after first solve)"),
    :is_linear => ("auto", "linear problem (avoid reassembly of nonlinear operators to check residual)"),
    :maxiterations => (10, "maximal number of nonlinear iterations/linear solves"),
    :method_linear => (UMFPACKFactorization(), "any solver or custom LinearSolveFunction compatible with LinearSolve.jl (default = UMFPACKFactorization())"),
    :plot => (false, "plot all solved unknowns with a (very rough but fast) unicode plot"),
    :precon_linear => (nothing, "function that computes preconditioner for method_linear in case an iterative solver is chosen"),
    :reltol => (1.0e-11, "reltol for linear solver (if iterative)"),
    :restrict_dofs => ([], "array of dofs for each unknown that should be solved (default: all dofs)"),
    :return_config => (false, "solver returns solver configuration (including A and b of last iteration)"),
    :show_config => (false, "show configuration at the beginning of solve"),
    :show_matrix => (false, "show system matrix after assembly"),
    :spy => (false, "show unicode spy plot of system matrix during solve"),
    :symmetrize => (false, "make system matrix symmetric (replace by (A+A')/2)"),
    :symmetrize_structure => (false, "make the system sparse matrix structurally symmetric (e.g. if [j,k] is also [k,j] must be set, all diagonal entries must be set)"),
    :target_residual => (1.0e-10, "stop if the absolute (nonlinear) residual is smaller than this number"),
    :time => (0.0, "current time to be used in all time-dependent operators"),
    :timeroutputs => (:full, "configures show of timeroutputs (choose between :hide, :full, :compact)"),
    :verbosity => (0, "verbosity level"),
)


function Base.show(io::IO, SC::SolverConfiguration)
    println(io, "\nSOLVER CONFIGURATION")
    println(io, "--------------------")
    println(io, "Problem:         ", SC.PD.name)
    println(io, "Unknowns:        ", join([u.name for u in SC.unknowns], ", "))
    println(io, "FE Spaces:       ", join([typeof(block.FES) for block in SC.sol], ", "))
    println(io, "Matrix size:     ", size(SC.A.entries))
    println(io, "Max iterations:  ", get(SC.parameters, :maxiterations, "N/A"))
    println(io, "Linear solver:   ", get(SC.parameters, :method_linear, "N/A"))
    println(io, "Preconditioner:  ", get(SC.parameters, :precon_linear, "N/A"))
    println(io, "Tolerance:       ", get(SC.parameters, :target_residual, "N/A"))
    println(io, "Constant matrix: ", get(SC.parameters, :constant_matrix, "N/A"))
    println(io, "Constant rhs:    ", get(SC.parameters, :constant_rhs, "N/A"))
    println(io, "Verbosity:       ", get(SC.parameters, :verbosity, "N/A"))
    println(io, "Initialized:     ", get(SC.parameters, :initialized, "N/A"))
    println(io, "--------------------")
    println(io, "Other parameters:")
    for (k, v) in SC.parameters
        if !(k in (:maxiterations, :method_linear, :precon_linear, :target_residual, :constant_matrix, :constant_rhs, :verbosity, :initialized))
            println(io, "  ", rpad(string(k), 16), " : ", v)
        end
    end
    return
end


"""
$(TYPEDSIGNATURES)

Construct a `SolverConfiguration` for a given problem and set of finite element spaces.

A `SolverConfiguration` bundles all data and options needed to assemble and solve a finite element system for a given `ProblemDescription`. It stores the system matrix, right-hand side, solution vector, solver parameters, and bookkeeping for unknowns and degrees of freedom.

# Arguments
- `Problem::ProblemDescription`: The problem to be solved, including operators, unknowns, and boundary conditions.
- `FES::Union{<:FESpace, Vector{<:FESpace}}`: The finite element space(s) for the unknowns. Can be a single space or a vector of spaces (one per unknown).
- `unknowns::Vector{Unknown}`: (optional) The unknowns to be solved for (default: `Problem.unknowns`).
- `init`: (optional) Initial `FEVector` for the solution. If provided, the finite element spaces are inferred from it.
- `kwargs...`: Additional keyword arguments to set solver parameters (see below).

# Keyword Arguments
$(_myprint(default_solver_kwargs()))

# Returns
- A `SolverConfiguration` instance, ready to be passed to the `solve` function.

# Notes
- If `init` is provided, the finite element spaces are inferred from it and the solution vector is initialized accordingly.
- The constructor checks that the number of unknowns matches the number of finite element spaces and that all unknowns are present in the problem description.
- The configuration includes storage for the system matrix, right-hand side, solution, temporary solution, residual, and solver statistics.

"""
function SolverConfiguration(Problem::ProblemDescription; init = nothing, unknowns = Problem.unknowns, kwargs...)
    ## try to guess FES from init
    if typeof(init) <: FEVector
        FES = [init[u].FES for u in unknowns]
    end
    return SolverConfiguration(Problem, unknowns, FES; kwargs...)
end

function SolverConfiguration(Problem::ProblemDescription, FES; unknowns = Problem.unknowns, kwargs...)
    return SolverConfiguration(Problem, unknowns, FES; kwargs...)
end

function SolverConfiguration(Problem::ProblemDescription, unknowns::Array{Unknown, 1}, FES, default_kwargs = default_solver_kwargs(); TvM = Float64, TiM = Int, bT = Float64, kwargs...)
    if typeof(FES) <: FESpace
        FES = [FES]
    end
    @assert length(unknowns) <= length(FES) "length of unknowns and FE spaces must coincide"
    ## check if unknowns are part of Problem description
    for u in unknowns
        @assert u in Problem.unknowns "unknown $u is not part of the given ProblemDescription"
    end
    parameters = Dict{Symbol, Any}(k => v[1] for (k, v) in default_kwargs)
    _update_params!(parameters, kwargs)
    ## compute offsets
    offsets = [0]
    for FE in FES
        push!(offsets, FE.ndofs + offsets[end])
    end

    ## storage for full system
    FES_active = FES[1:length(unknowns)]
    A = FEMatrix{TvM, TiM}(FES_active; tags = unknowns, npartitions = num_partitions(FES[1].xgrid))
    b = FEVector{bT}(FES_active; tags = unknowns)
    res = copy(b)

    ## initialize solution vector
    if parameters[:init] === nothing
        names = [u.name for u in unknowns]
        append!(names, ["N.N." for j in (length(unknowns) + 1):length(FES)])
        x = FEVector{bT}(FES; name = names, tags = unknowns)
        unknown_ids_in_sol = 1:length(unknowns)
    else
        x = parameters[:init]
        unknown_ids_in_sol = [findfirst(==(u), x.tags) for u in unknowns]
    end

    ## adjustments for using freedofs
    if haskey(parameters, :restrict_dofs)
        if length(parameters[:restrict_dofs]) > 0
            freedofs = Vector{Int}(parameters[:restrict_dofs][1])
            for j in 2:length(parameters[:restrict_dofs])
                parameters[:restrict_dofs][j] .+= FES[j - 1].ndofs
                append!(freedofs, parameters[:restrict_dofs][j])
            end
            x_temp = copy(b)
        else
            freedofs = []
            x_temp = x
        end
    else
        freedofs = []
        x_temp = x
    end

    return SolverConfiguration{typeof(A), typeof(b), typeof(x)}(Problem, A, b, x, x_temp, res, freedofs, default_statistics(TvM, TiM), nothing, unknown_ids_in_sol, unknowns, copy(unknowns), offsets, parameters)
end
