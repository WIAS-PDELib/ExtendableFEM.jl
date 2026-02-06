"""
    get_unknown_id(SC::SolverConfiguration, u::Unknown)

Returns the id of the unknown u in SC.
"""
function get_unknown_id(SC::SolverConfiguration, u::Unknown)
    return findfirst(==(u), SC.unknowns)
end
"""
    update_solution!(sol, x, unknowns, freedofs, damping)

Update the solution vector with the new values, including optional damping.
"""
function update_solution!(sol, x, unknowns, freedofs, damping)
    if length(freedofs) > 0
        sol.entries[freedofs] .= x
    else
        offset = 0
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

"""
    compute_nonlinear_residual!(residual, A, b, sol, unknowns, PD, SC, freedofs)

Compute the nonlinear residual for the current solution.
"""
function compute_nonlinear_residual!(residual, A, b, sol, unknowns, PD, SC, freedofs)
    residual.entries .= b.entries
    for j in 1:length(b), k in 1:length(b)
        addblock_matmul!(residual[j], A[j, k], sol[unknowns[k]], factor = -1.0)
    end

    # add Lagrange residuals
    for rs in PD.restrictions
        mul!(residual.entries, rs.parameters[:matrix], rs.parameters[:multiplier], -1.0, 1.0)
    end


    for op in PD.operators
        residual.entries[fixed_dofs(op)] .= 0
    end

    for u_off in SC.parameters[:inactive]
        j = get_unknown_id(SC, u_off)
        if j > 0
            fill!(residual[j], 0)
        end
    end

    nlres = length(freedofs) > 0 ? norm(residual.entries[freedofs]) : norm(residual.entries)
    restriction_residuals = [norm(rs.parameters[:matrix]' * sol.entries - rs.parameters[:rhs]) for rs in PD.restrictions]

    if length(PD.restrictions) > 0
        nlres = sqrt(nlres^2 + norm(restriction_residuals)^2)
    end

    for rs in PD.restrictions
        mul!(residual.entries, rs.parameters[:matrix], rs.parameters[:multiplier], 1.0, 1.0)
    end

    if SC.parameters[:verbosity] > 0
        if length(residual) > 1
            @info "sub-residuals = $(norms(residual))"
        end
        @info "nonlinear residuals of restrictions = $restriction_residuals"
    end

    return nlres
end

"""
    init_linear_solver!(SC, A, timer, method_linear, precon_linear)

Initialize the linear solver for the given system.
"""
function init_linear_solver!(SC, A, b, timer, method_linear, precon_linear)

    # TODO use the timer
    time_assembly = 0.0
    allocs_assembly = 0

    if SC.linsolver === nothing
        if SC.parameters[:verbosity] > 0
            @info ".... initializing linear solver ($(method_linear))\n"
        end
        @timeit timer "initialization" begin
            stats = @timed begin
                abstol = SC.parameters[:abstol]
                reltol = SC.parameters[:reltol]
                LP = LinearProblem(A, b)
                if precon_linear !== nothing
                    SC.linsolver = init(LP, method_linear; Pl = precon_linear(A), abstol, reltol)
                else
                    SC.linsolver = init(LP, method_linear; abstol, reltol)
                end
            end
            time_assembly += stats.time
            allocs_assembly += stats.bytes
        end
    end
    return time_assembly, allocs_assembly
end

"""
    assemble_system!(A, b, sol, PD, SC, timer; kwargs...)

Assemble the system matrix and right-hand side for the given problem.
"""
function assemble_system!(A, b, sol, PD, SC, timer; kwargs...)

    # TODO use the timer
    time_assembly = 0.0
    allocs_assembly = 0

    if !SC.parameters[:constant_rhs]
        fill!(b.entries, 0)
    end
    if !SC.parameters[:constant_matrix]
        fill!(A.entries.cscmatrix.nzval, 0)
    end

    # Assemble operators and restrictions
    if SC.parameters[:initialized]
        for op in PD.operators
            @timeit timer "$(op.parameters[:name])" begin
                stats = @timed assemble!(
                    A, b, sol, op, SC;
                    time = SC.parameters[:time],
                    assemble_matrix = !SC.parameters[:constant_matrix],
                    assemble_rhs = !SC.parameters[:constant_rhs],
                    kwargs...
                )
            end
            time_assembly += stats.time
            allocs_assembly += stats.bytes
        end
    else
        for op in PD.operators
            @timeit timer "$(op.parameters[:name]) (first)" begin
                stats = @timed assemble!(
                    A, b, sol, op, SC;
                    time = SC.parameters[:time],
                    kwargs...
                )
            end
            time_assembly += stats.time
            allocs_assembly += stats.bytes
        end
        ## assemble restrictions
        for restriction in PD.restrictions
            @timeit timer "$(restriction.parameters[:name])" assemble!(restriction, sol, SC; kwargs...)
        end
    end
    flush!(A.entries)

    # Apply penalties
    for op in PD.operators
        @timeit timer "$(op.parameters[:name]) (penalties)" begin
            stats = @timed apply_penalties!(
                A, b, sol, op, SC;
                assemble_matrix = !SC.parameters[:initialized] || !SC.parameters[:constant_matrix],
                assemble_rhs = !SC.parameters[:initialized] || !SC.parameters[:constant_rhs],
                kwargs...
            )
        end
        time_assembly += stats.time
        allocs_assembly += stats.bytes
    end
    flush!(A.entries)

    return time_assembly, allocs_assembly
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
    solve_linear_system!(A, b, sol, soltemp, residual, linsolve, unknowns, freedofs, damping, PD, SC, stats, is_linear, timer)

Solves the linear system and updates the solution vector. This includes:
- Setting up the system matrix and right-hand side
- Solving the linear system
- Computing the residual
- Updating the solution with optional damping
"""
function solve_linear_system!(A, b, sol, soltemp, residual, unknowns, freedofs, damping, PD, SC, stats, is_linear, timer, kwargs...)

    @timeit timer "linear solver" begin

        # does the linsolve object need a (new) matrix?
        linsolve_needs_matrix = !SC.parameters[:constant_matrix] || !SC.parameters[:initialized]

        ## start with the assembled matrix containing all assembled operators
        if linsolve_needs_matrix
            if length(freedofs) > 0
                A_unrestricted = A.entries.cscmatrix[freedofs, freedofs]
            else
                A_unrestricted = A.entries.cscmatrix
            end
        end

        # we solve for A Δx = r
        # and update x = sol - Δx
        if length(freedofs) > 0
            b_unrestricted = residual.entries[freedofs]
        else
            b_unrestricted = residual.entries
        end

        ## restrictions only involve the blocks coressponding to the unknowns and not necessarily the full sol.entries
        sol_freedofs = mortar([view(sol[u]) for u in unknowns])

        if length(PD.restrictions) == 0
            if linsolve_needs_matrix
                linsolve_A = A_unrestricted
            end
            linsolve_b = b_unrestricted
        else
            # add possible Lagrange restrictions
            restriction_matrices = [length(freedofs) > 0 ? view(restriction_matrix(re), freedofs, :) : restriction_matrix(re) for re in PD.restrictions ]
            restriction_rhss = deepcopy([length(freedofs) > 0 ? view(restriction_rhs(re), freedofs) : restriction_rhs(re) for re in PD.restrictions ])

            # block sizes for the block matrix
            block_sizes = [length(b_unrestricted), [ length(b) for b in restriction_rhss ]...]
            total_size = sum(block_sizes)

            # total number of additional LM dofs
            nLMs = @views sum(block_sizes[2:end])

            @timeit timer "LM restrictions (nLMs = $nLMs)" begin

                Tv = eltype(b_unrestricted)

                ## create block matrix
                if linsolve_needs_matrix
                    if SC.parameters[:compress_restrictions]
                        # combine all restriction matrices into one:
                        combined_transposed_restriction_matrix = vcat(restriction_matrices'...)
                        combined_restriction_rhs = vcat(restriction_rhss...)

                        # compute rank revealing QR decomposition (of the transposed matrix)
                        qr_result = qr(combined_transposed_restriction_matrix)

                        # extract components (from docs: Q*R = combined_transposed_restriction_matrix[prow, pcol])
                        (; Q, R, prow, pcol) = qr_result

                        # we need the inverse column permutation
                        ipcol = invperm(pcol)

                        # the new combined restriction block (add another transpose)
                        combined_restriction_matrix = R[:, ipcol]'

                        # the new combined restriction rhs
                        combined_restriction_rhs = Q'combined_restriction_rhs[prow]

                        # compress the  column space
                        qr_rank = rank(qr_result)
                        @assert norm(combined_restriction_rhs[(qr_rank + 1):end]) ≤ 1.0e-12 * norm(combined_restriction_rhs) "the rhs of the restriction is not in the image"

                        combined_restriction_matrix = combined_restriction_matrix[:, 1:qr_rank]
                        combined_restriction_rhs = combined_restriction_rhs[1:qr_rank]

                        # replace by single entries
                        restriction_matrices = [combined_restriction_matrix]
                        restriction_rhss = [combined_restriction_rhs]

                        # remove previously defined restrictions (with explicit copy of the rhs; it is otherwise modified later)
                        PD.restrictions = [CompressedRestriction(combined_restriction_matrix, deepcopy(combined_restriction_rhs))]

                        # update sizes
                        block_sizes = [length(b_unrestricted), length(combined_restriction_rhs)]
                        total_size = sum(block_sizes)
                    end

                    A_block = BlockMatrix(spzeros(Tv, total_size, total_size), block_sizes, block_sizes)
                    A_block[Block(1, 1)] = A_unrestricted
                end

                ## we need to add the (initial) solution to the rhs, since we work with the residual equation
                for (B, rhs) in zip(restriction_matrices, restriction_rhss)
                    rhs .-= B'sol_freedofs
                end

                b_block = BlockVector(zeros(Tv, total_size), block_sizes)
                b_block[Block(1)] = b_unrestricted

                for i in eachindex(restriction_matrices)
                    if linsolve_needs_matrix
                        A_block[Block(1, i + 1)] = restriction_matrices[i]
                        A_block[Block(i + 1, 1)] = transpose(restriction_matrices[i])
                    end
                    b_block[Block(i + 1)] = restriction_rhss[i]

                end

                # convert to dense vectors
                linsolve_b = Vector(b_block)

                if linsolve_needs_matrix

                    # linsolve.A = sparse(A_block) # convert to CSC Matrix is very slow https://github.com/JuliaArrays/BlockArrays.jl/issues/78
                    # do it manually:
                    A_flat = spzeros(size(A_block))

                    lasts_row = [0, axes(A_block)[1].lasts...] # add leading zero
                    lasts_col = [0, axes(A_block)[2].lasts...] # add leading zero

                    (n_row, n_col) = size(blocks(A_block))

                    for i in 1:n_row, j in 1:n_col
                        range_row = (lasts_row[i] + 1):lasts_row[i + 1]
                        range_col = (lasts_col[j] + 1):lasts_col[j + 1]

                        # write each block directly in the resulting matrix
                        A_flat[range_row, range_col] = A_block[Block(i, j)]
                    end
                    linsolve_A = A_flat
                end
            end
        end

        if !SC.parameters[:initialized]
            ## init solver if not done before
            @timeit timer "linear solver" begin
                init_linear_solver!(SC, linsolve_A, linsolve_b, timer, method_linear, precon_linear)
            end
            SC.parameters[:initialized] = true
        end

        linsolve = SC.linsolver

        # Solve linear system
        push!(stats[:matrix_nnz], nnz(linsolve.A))
        @timeit timer "solve! call" begin
            LinearSolve.solve!(linsolve)
        end
    end

    # Check linear residual
    @timeit timer "linear residual computation" begin

        # compute flat residual (reuse b_flat): residual_flat = A_flat * u_flat - b_flat
        residual_flat = linsolve.b
        mul!(residual_flat, linsolve.A, linsolve.u, 1.0, -1.0)

        for op in PD.operators
            for dof in fixed_dofs(op)
                # fix dofs only in first block
                if dof <= length(b_unrestricted)
                    residual_flat[dof] = 0
                end
            end
        end

        if length(freedofs) > 0
            residual.entries[freedofs] .= @views residual_flat[1:length(b_unrestricted)]
        else
            residual.entries .= @views residual_flat[1:length(b_unrestricted)]
        end

        if length(PD.restrictions) > 0
            # extract all residuals for the restriction blocks
            block_ends = cumsum(block_sizes)
            restriction_residuals = [norm(residual_flat[(block_ends[i] + 1):block_ends[i + 1]]) for i in 1:(length(block_sizes) - 1) ]
            push!(stats[:restriction_residuals], restriction_residuals)

            # store Lagrange multipliers
            for (i, rs) in enumerate(PD.restrictions)
                rs.parameters[:multiplier] = linsolve.u[(block_ends[i] + 1):block_ends[i + 1]]
            end
        end

        linres = norm(residual_flat)
        push!(stats[:linear_residuals], linres)
        if is_linear
            push!(stats[:nonlinear_residuals], linres)
        end
    end


    # Compute solution update
    @timeit timer "update solution" begin

        # extract the solution / dismiss the lagrange multipliers
        @views Δx = linsolve.u[1:length(b_unrestricted)]

        if length(freedofs) > 0
            x = sol.entries[freedofs] + Δx
        else
            x = zero(Δx)
            offset = 0
            for u in unknowns
                ndofs_u = length(view(sol[u]))
                x_range = (offset + 1):(offset + ndofs_u)
                x[x_range] .= view(sol[u]) .+ view(Δx, x_range)
                offset += ndofs_u
            end
        end
    end


    # Update solution
    @timeit timer "update solution" begin
        update_solution!(sol, x, unknowns, freedofs, damping)
    end

    return linres
end

"""
    initialize_coupled_solution(FES, init, unknowns, nPDs)

Initialize the solution vector and finite element spaces for coupled problems.
"""
function initialize_coupled_solution(FES, init, unknowns, nPDs)
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
    return sol, FES
end

"""
    solve_coupled_system!(A, b, sol, residual, linsolve, unknowns, damping, PD, SC)

Solves the coupled system for a single subproblem and updates the solution.
"""
function solve_coupled_system!(A, b, sol, residual, linsolve, unknowns, damping, PD, SC)
    if !SC.parameters[:constant_matrix] || !SC.parameters[:initialized]
        linsolve.A = A.entries.cscmatrix
    end

    # we solve for A Δx = r
    # and update x = sol - Δx
    linsolve.b = residual.entries
    SC.parameters[:initialized] = true

    ## solve
    result = LinearSolve.solve!(linsolve)
    Δx = result.u

    # x = sol.entries - Δx ... in the entry ranges of the present unknowns
    x = zero(Δx)
    offset = 0
    for u in unknowns
        ndofs_u = length(view(sol[u]))
        x_range = (offset + 1):(offset + ndofs_u)
        x[x_range] .= view(sol[u]) .- view(Δx, x_range)
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

    linres = norm(residual.entries)

    offset = 0
    for u in unknowns
        ndofs_u = length(view(sol[u]))
        if damping > 0
            view(sol[u]) .= damping * view(sol[u]) + (1 - damping) * view(x, (offset + 1):(offset + ndofs_u))
        else
            view(sol[u]) .= view(x, (offset + 1):(offset + ndofs_u))
        end
        offset += ndofs_u
    end

    return linres
end

"""
    check_problem_linearity!(PDs, SCs, unknowns, is_linear)

Check the linearity of each subproblem and set appropriate flags.
"""
function check_problem_linearity!(PDs, SCs, unknowns)
    nPDs = length(PDs)
    is_linear = zeros(Bool, nPDs)
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
        elseif SCs[j].parameters[:is_linear] == true
            is_linear[j] = true
        end
        if is_linear[j] && nonlinear[j]
            @warn "problem $(PD.name) seems nonlinear, but user set is_linear = true (results may be wrong)!!"
        end
    end
    return is_linear
end

"""
    init_solver_config(PD::ProblemDescription, FES::Union{<:FESpace, Vector{<:FESpace}}, SC, unknowns, kwargs)

Initialize and configure the solver based on the problem description and finite element spaces.
"""
function init_solver_config(PD::ProblemDescription, FES::Union{<:FESpace, Vector{<:FESpace}}, SC, unknowns, kwargs)
    FES_array = typeof(FES) <: FESpace ? [FES] : FES

    if typeof(SC) <: SolverConfiguration
        _update_params!(SC.parameters, kwargs)
        SC.parameters[:verbosity] > 0 && @info ".... reusing given solver configuration\n"
    else
        SC = SolverConfiguration(PD, unknowns, FES_array; kwargs...)
        SC.parameters[:verbosity] > 0 && @info ".... init solver configuration\n"
    end
    return SC, FES_array
end


function print_convergence_result(SC, is_linear, linres, nlres, nltol, j, maxits)

    stop = false
    if nlres < nltol
        if SC.parameters[:verbosity] > -1
            @printf " END\t"
            @printf "%.3e\t" nlres
            @printf "converged\n"
        end
        stop = true
    elseif isnan(nlres)
        if SC.parameters[:verbosity] > -1
            @printf " END\t"
            @printf "%.3e\t" nlres
            @printf "not converged\n"
        end
        stop = true
    elseif (j == maxits + 1) && !(is_linear)
        if SC.parameters[:verbosity] > -1
            @printf " END\t"
            @printf "\t\t%.3e\t" linres
            @printf "maxiterations reached\n"
        end
        stop = true
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
    return stop
end

function CommonSolve.solve(PD::ProblemDescription, FES::Union{<:FESpace, Vector{<:FESpace}}, SC = nothing; unknowns = PD.unknowns, kwargs...)
    SC, FES = init_solver_config(PD, FES, SC, unknowns, kwargs)

    ## load TimerOutputs
    timer = timeroutputs(SC)
    if timer == Any[]
        timer = TimerOutput()
    end

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
    is_linear = check_problem_linearity!([PD], [SC], unknowns)[1]
    if is_linear
        maxits = 0
    end

    if SC.parameters[:verbosity] > -1
        @printf " #IT\t------- RESIDUALS -------\n"
        @printf "   \tNONLINEAR\tLINEAR\n"
    end
    nlres = 1.1e30
    linres = 1.1e30

    for j in 1:(maxits + 1)
        if is_linear && j == 2
            nlres = linres
        else
            @timeit timer "assembly" begin
                assemble_system!(A, b, sol, PD, SC, timer; kwargs...)
            end

            ## show spy
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

            ## compute nonlinear residual
            @timeit timer "assembly" @timeit timer "residual vector" begin
                nlres = compute_nonlinear_residual!(residual, A, b, sol, unknowns, PD, SC, freedofs)
            end
        end
        if !is_linear
            push!(stats[:nonlinear_residuals], nlres)
        end

        stop = print_convergence_result(SC, is_linear, linres, nlres, nltol, j, maxits)
        if stop
            break
        end

        linres = solve_linear_system!(
            A, b, sol, soltemp, residual, unknowns,
            freedofs, damping, PD, SC, stats, is_linear, timer
        )

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

    # Initialize solution vector and FE spaces
    sol, FES = initialize_coupled_solution(FES, init, unknowns, nPDs)

    @info "SOLVING iteratively $([PD.name for PD in PDs])
			unknowns = $([[uj.name for uj in u] for u in unknowns])"
    #      fetypes = $(["$(get_FEType(FES[j]))" for j = 1 : length(unknowns)])
    #      ndofs = $([FES[j].ndofs for j = 1 : length(unknowns)])

    As = [SC.A for SC in SCs]
    bs = [SC.b for SC in SCs]
    residuals = [SC.res for SC in SCs]

    # Check linearity of each subproblem
    is_linear = check_problem_linearity!(PDs, SCs, unknowns)
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
                    # Assemble system and update timing/allocation info
                    assembly_time, assembly_allocs = assemble_system!(A, b, sol, PD, SC, TimerOutput(); kwargs...)
                    time_assembly += assembly_time
                    allocs_assembly += assembly_allocs

                    # Initialize linear solver if needed
                    solve_init_time, solve_init_allocs = init_linear_solver!(
                        SC,
                        A,
                        TimerOutput(),
                        SC.parameters[:method_linear],
                        SC.parameters[:precon_linear]
                    )
                    time_solve_init += solve_init_time
                    allocs_solve_init += solve_init_allocs

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
                        linres = solve_coupled_system!(A, b, sol, residual, SC.linsolver, unknowns[p], damping, PD, SC)
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
