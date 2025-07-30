"""
````
get_unknown_id(SC::SolverConfiguration, u::Unknown)
````

returns the id of the unknown u in SC

"""
function get_unknown_id(SC::SolverConfiguration, u::Unknown)
    return findfirst(==(u), SC.unknowns)
end


"""
$(TYPEDSIGNATURES)

Solve a PDE problem described by a `ProblemDescription` using the provided finite element spaces and solver configuration.

# Arguments
- `PD::ProblemDescription`: The problem description, including operators, unknowns, and boundary conditions.
- `FES::Union{<:FESpace, Vector{<:FESpace}, Dict{<:Unknown}}`: The finite element space(s) for discretizing the unknowns. Can be a single space, a vector of spaces (one per unknown), or a dictionary mapping unknowns to spaces.
- `SC`: (optional) Solver configuration. If not provided, a new configuration is created.
- `unknowns`: (optional) Vector of unknowns to solve for (default: `PD.unknowns`).
- `kwargs...`: Additional keyword arguments for solver configuration and algorithmic options (see below).

# Keyword Arguments
$(_myprint(default_solver_kwargs()))

# Returns
- The solution as an `FEVector` (or a tuple `(sol, SC)` if `return_config=true`).

# Notes
- The function automatically detects and handles linear and nonlinear problems depending on the operator definitions and runs a fixed-point iteration if necessary.
- See the `iterate_until_stationarity` function for non-monolithic problems.
- This method extends the `CommonSolve.solve` interface, where `ProblemDescription` acts as the problem type and `FES` as the solver type.

"""
function CommonSolve.solve(PD::ProblemDescription, FES::Dict{<:Unknown}, SC = nothing; unknowns = PD.unknowns, kwargs...)
    return solve(PD, [FES[u] for u in unknowns], SC; unknowns = unknowns, kwargs...)
end

function CommonSolve.solve(PD::ProblemDescription, SC = nothing; init = nothing, unknowns = PD.unknowns, kwargs...)
    if init === nothing
        @error "need to know initial FEVector or finite element spaces for unknowns of problem"
    end
    return solve(PD, [init[u].FES for u in unknowns], SC; init = init, unknowns = unknowns, kwargs...)
end

function symmetrize_structure!(A::ExtendableSparseMatrix{Tv, Ti}; diagval = 1.0e-16) where {Tv, Ti}
    cscmat = A.cscmatrix
    rows::Array{Ti, 1} = rowvals(cscmat)
    nzvals::Array{Tv, 1} = cscmat.nzval
    for col in 1:size(cscmat, 2)
        for r in nzrange(cscmat, col)
            A[col, rows[r]] += 1.0e-16 * nzvals[r]
        end
        if A[col, col] == 0 && diagval !== 0
            A[col, col] = diagval
        end
    end
    return flush!(A)
end

"""
    SolverDetails

    Module containing internal solve-subroutines.
"""
module SolverDetails

    function check_nonlinear(PD, unknowns)
        for op in PD.operators
            nl_dependencies = depends_nonlinearly_on(op)
            for u in unknowns
                if u in nl_dependencies
                    return true
                end
            end
        end
        return false
    end

    ## assemble operators
    function assembly!(A, b, sol, SC, PD, kwargs...)
        if !SC.parameters[:constant_rhs]
            fill!(b.entries, 0)
        end
        if !SC.parameters[:constant_matrix]
            fill!(A.entries.cscmatrix.nzval, 0)
        end
        if SC.parameters[:initialized]
            for op in PD.operators
                @timeit timer "$(op.parameters[:name])" assemble!(A, b, sol, op, SC; time = SC.parameters[:time], assemble_matrix = !SC.parameters[:constant_matrix], assemble_rhs = !SC.parameters[:constant_rhs], kwargs...)
            end
        else
            for op in PD.operators
                @timeit timer "$(op.parameters[:name]) (first)" assemble!(A, b, sol, op, SC; time = SC.parameters[:time], kwargs...)
            end
        end
        flush!(A.entries)

        ## penalize fixed dofs
        for op in PD.operators
            @timeit timer "$(op.parameters[:name]) (penalties)" apply_penalties!(A, b, sol, op, SC; assemble_matrix = !SC.parameters[:initialized] || !SC.parameters[:constant_matrix], assemble_rhs = !SC.parameters[:initialized] || !SC.parameters[:constant_rhs], kwargs...)
        end
        return flush!(A.entries)
    end

    function init_linsolve!(SC, linsolve, method_linear, abstol, reltol)
        if linsolve === nothing
            if SC.parameters[:verbosity] > 0
                @info ".... initializing linear solver ($(method_linear))\n"
            end
            @timeit timer "initialization" begin
                abstol = SC.parameters[:abstol]
                reltol = SC.parameters[:reltol]
                if precon_linear !== nothing
                    linsolve = init(SC.LP, method_linear; Pl = precon_linear(A.entries.cscmatrix), abstol, reltol)
                else
                    linsolve = init(SC.LP, method_linear; abstol, reltol)
                end
                SC.linsolver = linsolve
            end
        end
        return linsolve
    end

    function compute_nonlinear_residual!(residual, A, b, sol, PD, SC, freedofs)
        fill!(residual.entries, 0)
        for j in 1:length(b), k in 1:length(b)
            addblock_matmul!(residual[j], A[j, k], sol[unknowns[k]])
        end
        residual.entries .-= b.entries
        for op in PD.operators
            residual.entries[fixed_dofs(op)] .= 0
        end
        for u_off in SC.parameters[:inactive]
            j = get_unknown_id(SC, u_off)
            if j > 0
                fill!(residual[j], 0)
            end
        end
        if length(freedofs) > 0
            nlres = norm(residual.entries[freedofs])
        else
            nlres = norm(residual.entries)
        end
        if SC.parameters[:verbosity] > 0 && length(residual) > 1
            @info "sub-residuals = $(norms(residual))"
        end
        return nothing
    end

    function show_spy!(A, SC)
        if SC.parameters[:symmetrize]
            A.entries.cscmatrix = (A.entries.cscmatrix + A.entries.cscmatrix') / 2
        elseif SC.parameters[:symmetrize_structure]
            symmetrize_structure!(A.entries)
        end
        if SC.parameters[:show_matrix]
            @show A
        elseif SC.parameters[:spy]
            @info ".... spy plot of system matrix:\n$(A.entries.cscmatrix))"
        end
        if SC.parameters[:check_matrix]
            @info ".... ||A - A'|| = $(norm(A.entries.cscmatrix - A.entries.cscmatrix', Inf))"
            @info "....  isposdef  = $(isposdef(A.entries.cscmatrix))"
        end
        return nothing
    end

    function print_stats!(stats, nlres, nltol, SC, j, maxits, is_linear, linres)
        if !is_linear
            push!(stats[:nonlinear_residuals], nlres)
        end

        if nlres < nltol
            if SC.parameters[:verbosity] > -1
                @printf " END\t"
                @printf "%.3e\t" nlres
                @printf "converged\n"
            end
            break
        elseif isnan(nlres)
            if SC.parameters[:verbosity] > -1
                @printf " END\t"
                @printf "%.3e\t" nlres
                @printf "not converged\n"
            end
            break
        elseif (j == maxits + 1) && !(is_linear)
            if SC.parameters[:verbosity] > -1
                @printf " END\t"
                @printf "\t\t%.3e\t" linres
                @printf "maxiterations reached\n"
            end
            break
        else
            if SC.parameters[:verbosity] > -1
                if is_linear
                    @printf " END\t"
                else
                    @printf "%4d\t" j
                end
                if !(is_linear)
                    @printf "%.3e\t" nlres
                else
                    @printf "---------\t"
                end
            end
        end

        return nothing
    end

    function prepare_solution!(sol, Δx, freedofs, unknowns)
        if length(freedofs) > 0
            x = sol.entries[freedofs] - Δx.u
        else
            x = zero(Δx)
            offset = 0
            for u in unknowns
                ndofs_u = length(view(sol[u]))
                x_range = (offset + 1):(offset + ndofs_u)
                x[x_range] .= view(sol[u]) .- view(Δx, x_range)
                offset += ndofs_u
            end
        end
        return nothing
    end

    function compute_linear_residual!(residual, A, b, x, soltemp, PD)
        if length(freedofs) > 0
            soltemp.entries[freedofs] .= x
            residual.entries .= A.entries.cscmatrix * soltemp.entries
        else
            residual.entries .= A.entries.cscmatrix * x
        end
        residual.entries .-= b.entries
        for op in PD.operators
            for dof in fixed_dofs(op)
                if dof <= length(residual.entries)
                    residual.entries[dof] = 0
                end
            end
        end
        return nothing
    end

    function update_solution!(sol, x, freedofs, unknowns)
        offset = 0
        if length(freedofs) > 0
            sol.entries[freedofs] .= x
        else
            for u in unknowns
                ndofs_u = length(view(sol[u]))
                if damping > 0
                    view(sol[u]) .= damping * view(sol[u]) + (1 - damping) * view(x, (offset + 1):(offset + ndofs_u))
                else
                    view(sol[u]) .= view(x, (offset + 1):(offset + ndofs_u))
                end
                offset += ndofs_u
            end
        end
        return nothing
    end
end


function CommonSolve.solve(PD::ProblemDescription, FES::Union{<:FESpace, Vector{<:FESpace}}, SC = nothing; unknowns = PD.unknowns, kwargs...)

    if typeof(FES) <: FESpace
        FES = [FES]
    end

    if typeof(SC) <: SolverConfiguration
        _update_params!(SC.parameters, kwargs)
        SC.parameters[:verbosity] > 0 && @info ".... reusing given solver configuration\n"
    else
        SC = SolverConfiguration(PD, unknowns, FES; kwargs...)
        SC.parameters[:verbosity] > 0 && @info ".... init solver configuration\n"
    end

    ## load TimerOutputs
    timer = timeroutputs(SC)
    if timer == Any[]
        timer = TimerOutput()
    end

    ## prepare some data
    A = SC.A
    b = SC.b
    sol = SC.sol
    soltemp = SC.tempsol
    residual = SC.res
    method_linear = SC.parameters[:method_linear]
    precon_linear = SC.parameters[:precon_linear]
    stats = SC.statistics
    for (key, value) in stats
        stats[key] = []
    end

    ## unpack solver parameters
    maxits = SC.parameters[:maxiterations]
    @assert maxits > -1
    nltol = SC.parameters[:target_residual]
    is_linear = SC.parameters[:is_linear]
    damping = SC.parameters[:damping]
    freedofs = SC.freedofs

    if SC.parameters[:verbosity] > -1
        if length(freedofs) > 0
            @info "SOLVING $(PD.name) @ time = $(SC.parameters[:time])
				unknowns = $([u.name for u in unknowns])
				fetypes = $(["$(get_FEType(FES[j]))" for j in 1:length(unknowns)])
				ndofs = $([FES[j].ndofs for j in 1:length(unknowns)]) (restricted to $(length.(SC.parameters[:restrict_dofs])))"
        else
            @info "SOLVING $(PD.name) @ time = $(SC.parameters[:time])
				unknowns = $([u.name for u in unknowns])
				fetypes = $(["$(get_FEType(FES[j]))" for j in 1:length(unknowns)])
				ndofs = $([FES[j].ndofs for j in 1:length(unknowns)])"
        end
    end

    if SC.parameters[:verbosity] > 0 || SC.parameters[:show_config]
        @info "\n$(SC)"
    end

    ## check if problem is (non)linear
    nonlinear = SolverDetails.check_nonlinear(PD, unknowns)

    SC.parameters[:verbosity] > 0 && @info " nonlinear = $(nonlinear ? "true" : "false")\n"
    is_linear == "auto" && (is_linear = !nonlinear)

    if is_linear && nonlinear
        @warn "problem seems nonlinear, but user set is_linear = true (results may be wrong)!!"
    end

    is_linear && (maxits = 0)

    if SC.parameters[:verbosity] > -1
        @printf " #IT\t------- RESIDUALS -------\n"
        @printf "   \tNONLINEAR\tLINEAR\n"
    end

    nlres = 1.1e30
    linres = 1.1e30
    linsolve = SC.linsolver

    for j in 1:(maxits + 1)
        if is_linear && j == 2
            nlres = linres
        else
            @timeit timer "assembly" SolverDetails.assembly!(A, b, sol, SC, PD, kwargs...)

            ## show spy
            SolverDetails.show_spy!(A, SC)

            ## init solver
            @timeit timer "linear solver" begin
                linsolve = init_linsolve!(SC, linsolve, method_linear, abstol, reltol)
            end

            ## compute nonlinear residual
            @timeit timer "assembly" @timeit timer "residual vector" begin
                compute_nonlinear_residual!(residual, A, b, sol, PD, SC, freedofs)
            end
        end


        SolverDetails.print_stats!(stats, nlres, nltol, SC, j, maxits, is_linear, linres)

        @timeit timer "linear solver" begin

            if !SC.parameters[:constant_matrix] || !SC.parameters[:initialized]
                if length(freedofs) > 0
                    linsolve.A = A.entries.cscmatrix[freedofs, freedofs]
                else
                    linsolve.A = A.entries.cscmatrix
                end
            end

            # we solve for A Δx = r
            # and update x = sol - Δx
            if length(freedofs) > 0
                linsolve.b = residual.entries[freedofs]
            else
                linsolve.b = residual.entries
            end
            SC.parameters[:initialized] = true

            ## solve
            push!(stats[:matrix_nnz], nnz(linsolve.A))
            @timeit timer "solve! call" Δx = LinearSolve.solve!(linsolve)

            # x = sol.entries - Δx.u for free dofs or partial solutions
            @timeit timer "update solution" begin
                SolverDetails.prepare_solution!(sol, Δx, freedofs, unknowns)
            end

            ## check linear residual with full matrix
            @timeit timer "linear residual computation" begin
                SolverDetails.compute_linear_residual!(residual, A, b, x, soltemp, PD)
            end

            linres = norm(residual.entries)
            push!(stats[:linear_residuals], linres)
            if is_linear
                push!(stats[:nonlinear_residuals], linres)
            end

            ## update solution (incl. damping etc.)
            @timeit timer "update solution" begin
                SolverDetails.update_solution!(sol, x, freedofs, unknowns)
            end
        end


        if SC.parameters[:verbosity] > -1
            if is_linear
                @printf "%.3e\t" linres
                @printf "finished\n"
            else
                @printf "%.3e\n" linres
            end
        end
    end

    if SC.parameters[:plot]
        for u in unknowns
            println(stdout, unicode_scalarplot(sol[u]; title = u.name, kwargs...))
        end
    end

    # Print the timings in the default way
    stats[:timeroutputs] = timer
    if SC.parameters[:timeroutputs] in [:full, :compact]
        # TimerOutputs.complement!(to)
        @printf "\n"
        print_timer(timer, compact = SC.parameters[:timeroutputs] == :compact, sortby = :firstexec)
    end

    if SC.parameters[:return_config]
        return sol, SC
    else
        return sol
    end
end


"""
$(TYPEDSIGNATURES)

Iteratively solve coupled problems (i.e. an array of SolverConfigurations/ProblemDescriptions with intersecting unknowns) by iterating between them. Iteration continues until the residuals of all subproblems fall below their specified tolerances, or until `maxsteps` is reached.

# Arguments
- `SCs::Vector{<:SolverConfiguration}`: Array of solver configurations, one for each subproblem. Each must contain its own `ProblemDescription`.
- `FES`: (optional) Nested vector of finite element spaces for all subproblems and unknowns. If not provided, must be inferred from `init`.
- `maxsteps::Int`: Maximum number of outer iterations (default: 1000).
- `energy_integrator`: (optional) Function or object to evaluate an energy or error after each iteration (for diagnostics).
- `init`: (optional) Initial `FEVector` containing all unknowns for all subproblems. Required if `FES` is not provided.
- `unknowns`: (optional) List of unknowns for each subproblem (default: `[SC.PD.unknowns for SC in SCs]`).
- `kwargs...`: Additional keyword arguments passed to each subproblem's solver.

# Returns
- `sol`: The final solution vector containing all unknowns.
- `it`: The number of outer iterations performed.

# Notes
- Each subproblem is solved in sequence within each outer iteration, using its own configuration and options.´
- If `energy_integrator` is provided, its value is printed after each iteration for monitoring convergence.
- This function is intended for non-monolithic (partitioned) solution strategies; for monolithic problems with a single ProblemDescription, see `solve`.


"""
function iterate_until_stationarity(
        SCs::Array{<:SolverConfiguration, 1},
        FES = nothing;
        maxsteps = 1000,
        energy_integrator = nothing,
        init = nothing,
        unknowns = [SC.PD.unknowns for SC in SCs],
        kwargs...
    )

    PDs::Array{ProblemDescription, 1} = [SC.PD for SC in SCs]
    nPDs = length(PDs)

    ## find FESpaces and generate solution vector
    if FES === nothing
        @assert init !== nothing "need init vector or FES (as a Vector{Vector{<:FESpace}})"
        @info ".... taking FESpaces from init vector \n"
        all_unknowns = init.tags
        for p in 1:nPDs, u in unknowns[p]
            @assert u in all_unknowns "did not found unknown $u in init vector (tags missing?)"
        end
        FES = [[init[u].FES for u in unknowns[j]] for j in 1:nPDs]
        sol = copy(init)
        sol.tags .= init.tags
    else
        all_unknowns = []
        for p in 1:nPDs, u in unknowns[p]
            if !(u in all_unknowns)
                push!(u, all_unknowns)
            end
        end
        sol = FEVector(FES; tags = all_unknowns)
    end

    @info "SOLVING iteratively $([PD.name for PD in PDs])
			unknowns = $([[uj.name for uj in u] for u in unknowns])"
    #      fetypes = $(["$(get_FEType(FES[j]))" for j = 1 : length(unknowns)])
    #      ndofs = $([FES[j].ndofs for j = 1 : length(unknowns)])

    As = [SC.A for SC in SCs]
    bs = [SC.b for SC in SCs]
    residuals = [SC.res for SC in SCs]

    ## unpack solver parameters
    is_linear = zeros(Bool, nPDs)

    ## check if problems are (non)linear
    nonlinear = zeros(Bool, nPDs)
    for (j, PD) in enumerate(PDs)
        for op in PD.operators
            nl_dependencies = depends_nonlinearly_on(op)
            for u in unknowns
                if u in nl_dependencies
                    nonlinear[j] = true
                    break
                end
            end
        end
        if SCs[j].parameters[:verbosity] > 0
            @info "nonlinear = $(nonlinear[j] ? "true" : "false")\n"
        end
        if SCs[j].parameters[:is_linear] == "auto"
            is_linear[j] = !nonlinear[j]
        end
        if is_linear[j] && nonlinear[j]
            @warn "problem $(PD.name) seems nonlinear, but user set is_linear = true (results may be wrong)!!"
        end
    end
    maxits = [is_linear[j] ? 1 : maxits[j] for j in 1:nPDs]

    alloc_factor = 1024^2

    time_final = 0
    allocs_final = 0
    nlres = 1.1e30
    linres = 1.1e30
    converged = zeros(Bool, nPDs)
    it::Int = 0
    while (it < maxsteps) && (any(converged .== false))
        it += 1
        @printf "%5d\t" it
        copyto!(init.entries, sol.entries)
        allocs_assembly = 0
        time_assembly = 0
        time_total = 0
        time_solve_init = 0
        allocs_solve_init = 0
        for p in 1:nPDs
            b = bs[p]
            A = As[p]
            PD = PDs[p]
            SC = SCs[p]
            residual = residuals[p]
            maxits = SC.parameters[:maxiterations]
            nltol = SC.parameters[:target_residual]
            damping = SC.parameters[:damping]
            for j in 1:1
                time_total += @elapsed begin

                    ## assemble operators
                    if !SC.parameters[:constant_rhs]
                        fill!(b.entries, 0)
                    end
                    if !SC.parameters[:constant_matrix]
                        fill!(A.entries.cscmatrix.nzval, 0)
                    end
                    if SC.parameters[:initialized]
                        time_assembly += @elapsed for op in PD.operators
                            allocs_assembly += @allocated assemble!(A, b, sol, op, SC; time = SC.parameters[:time], assemble_matrix = !SC.parameters[:constant_matrix], assemble_rhs = !SC.parameters[:constant_rhs], kwargs...)
                        end
                    else
                        time_assembly += @elapsed for op in PD.operators
                            allocs_assembly += @allocated assemble!(A, b, sol, op, SC; time = SC.parameters[:time], kwargs...)
                        end
                    end
                    flush!(A.entries)

                    ## penalize fixed dofs
                    time_assembly += @elapsed for op in PD.operators
                        allocs_assembly += @allocated apply_penalties!(A, b, sol, op, SC; kwargs...)
                    end
                    flush!(A.entries)

                    if SC.parameters[:verbosity] > 0
                        @printf " assembly time | allocs = %.2f s | %.2f MiB\n" time allocs / alloc_factor
                    end

                    ## show spy
                    if SC.parameters[:show_matrix]
                        @show A
                    elseif SC.parameters[:spy]
                        @info ".... spy plot of system matrix:\n$(UnicodePlots.spy(sparse(A.entries.cscmatrix)))"
                    end

                    ## init solver
                    linsolve = SC.linsolver
                    if linsolve === nothing
                        if SC.parameters[:verbosity] > 0
                            @info ".... initializing linear solver ($(method_linear))\n"
                        end
                        time_solve_init += @elapsed begin
                            allocs_solve_init += @allocated begin
                                method_linear = SC.parameters[:method_linear]
                                precon_linear = SC.parameters[:precon_linear]
                                abstol = SC.parameters[:abstol]
                                reltol = SC.parameters[:reltol]
                                LP = SC.LP
                                if precon_linear !== nothing
                                    linsolve = LinearSolve.init(LP, method_linear; Pl = precon_linear(linsolve.A), abstol = abstol, reltol = reltol)
                                else
                                    linsolve = LinearSolve.init(LP, method_linear; abstol = abstol, reltol = reltol)
                                end
                                SC.linsolver = linsolve
                            end
                        end
                    end

                    ## compute nonlinear residual
                    fill!(residual.entries, 0)
                    for j in 1:length(b), k in 1:length(b)
                        addblock_matmul!(residual[j], A[j, k], sol[unknowns[p][k]])
                    end
                    residual.entries .-= b.entries
                    #res = A.entries * sol.entries - b.entries
                    for op in PD.operators
                        residual.entries[fixed_dofs(op)] .= 0
                    end
                    for u_off in SC.parameters[:inactive]
                        j = get_unknown_id(SC, u_off)
                        if j > 0
                            fill!(residual[j], 0)
                        end
                    end
                    nlres = norm(residual.entries)
                    @printf "\tres[%d] = %.2e" p nlres
                end
                time_final += time_assembly + time_solve_init
                allocs_final += allocs_assembly + allocs_solve_init

                if nlres < nltol
                    converged[p] = true
                else
                    converged[p] = false
                end

                time_solve = @elapsed begin
                    allocs_solve = @allocated begin
                        if !SC.parameters[:constant_matrix] || !SC.parameters[:initialized]
                            linsolve.A = A.entries.cscmatrix
                        end

                        # we solve for A Δx = r
                        # and update x = sol - Δx
                        linsolve.b = residual.entries
                        SC.parameters[:initialized] = true


                        ## solve
                        Δx = LinearSolve.solve!(linsolve)

                        # x = sol.entries - Δx.u ... in the entry ranges of the present unknowns
                        x = zero(Δx.u)
                        offset = 0
                        for u in unknowns[p]
                            ndofs_u = length(view(sol[u]))
                            x_range = (offset + 1):(offset + ndofs_u)
                            x[x_range] .= view(sol[u]) .- view(Δx.u, x_range)
                            offset += ndofs_u
                        end

                        fill!(residual.entries, 0)
                        mul!(residual.entries, A.entries.cscmatrix, x)
                        residual.entries .-= b.entries
                        for op in PD.operators
                            for dof in fixed_dofs(op)
                                if dof <= length(residual.entries)
                                    residual.entries[dof] = 0
                                end
                            end
                        end
                        #@info residual.entries, norms(residual)
                        linres = norm(residual.entries)
                        offset = 0
                        for u in unknowns[p]
                            ndofs_u = length(view(sol[u]))
                            if damping > 0
                                view(sol[u]) .= damping * view(sol[u]) + (1 - damping) * view(x, (offset + 1):(offset + ndofs_u))
                            else
                                view(sol[u]) .= view(x, (offset + 1):(offset + ndofs_u))
                            end
                            offset += ndofs_u
                        end
                    end
                end
                time_total += time_solve
                time_final += time_solve
                allocs_final += allocs_solve
                time_solve += time_solve_init
                allocs_solve += allocs_solve_init
                if SC.parameters[:verbosity] > -1
                    @printf " (%.3e)" linres
                end
            end # nonlinear iterations subproblem
        end

        if energy_integrator !== nothing
            error = evaluate(energy_integrator, sol)
            @printf "   energy = %.3e" sum([sum(view(error, j, :)) for j in 1:size(error, 1)])
        end
        @printf "\n"
    end

    return sol, it
end
