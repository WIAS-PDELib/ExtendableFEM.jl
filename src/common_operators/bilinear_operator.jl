mutable struct BilinearOperatorFromMatrix{UT <: Union{Unknown, Integer}, MT} <: AbstractOperator
    u_test::Array{UT, 1}
    u_ansatz::Array{UT, 1}
    u_args::Array{UT, 1}
    A::MT
    parameters::Dict{Symbol, Any}
end

# informs solver when operator needs reassembly
function depends_nonlinearly_on(O::BilinearOperatorFromMatrix)
    return O.u_args
end

# informs solver in which blocks the operator assembles to
function dependencies_when_linearized(O::BilinearOperatorFromMatrix)
    return [O.u_ansatz, O.u_test]
end

# informs solver when operator needs reassembly in a time dependent setting
function is_timedependent(O::BilinearOperatorFromMatrix)
    return O.parameters[:time_dependent]
end

function Base.show(io::IO, O::BilinearOperatorFromMatrix)
    dependencies = dependencies_when_linearized(O)
    print(io, "$(O.parameters[:name])($([ansatz_function(dependencies[1][j]) for j in 1:length(dependencies[1])]), $([test_function(dependencies[2][j]) for j in 1:length(dependencies[2])]))")
    return nothing
end


mutable struct BilinearOperator{Tv <: Real, UT <: Union{Unknown, Integer}, KFT, MT} <: AbstractOperator
    u_test::Array{UT, 1}
    ops_test::Array{DataType, 1}
    u_ansatz::Array{UT, 1}
    ops_ansatz::Array{DataType, 1}
    u_args::Array{UT, 1}
    ops_args::Array{DataType, 1}
    kernel::KFT
    BE_test_vals::Array{Array{Array{Tv, 3}, 1}}
    BE_ansatz_vals::Array{Array{Array{Tv, 3}, 1}}
    BE_args_vals::Array{Array{Array{Tv, 3}, 1}}
    FES_test::Any             #::Array{FESpace,1}
    FES_ansatz::Any           #::Array{FESpace,1}
    FES_args::Any             #::Array{FESpace,1}
    BE_test::Any              #::Union{Nothing, Array{FEEvaluator,1}}
    BE_ansatz::Any            #::Union{Nothing, Array{FEEvaluator,1}}
    BE_args::Any              #::Union{Nothing, Array{FEEvaluator,1}}
    QP_infos::Any             #::Array{QPInfosT,1}
    L2G::Any
    QF::Any
    assembler::Any
    storage::MT
    parameters::Dict{Symbol, Any}
end

default_blfop_kwargs() = Dict{Symbol, Tuple{Any, String}}(
    :entities => (ON_CELLS, "assemble operator on these grid entities (default = ON_CELLS)"),
    :name => ("BilinearOperator", "name for operator used in printouts"),
    :transposed_copy => (0, "assemble a transposed copy of that operator into the transposed matrix block(s), 0 = no, 1 = symmetric, -1 = skew-symmetric"),
    :factor => (1, "factor that should be multiplied during assembly"),
    :lump => (0, "diagonal lumping (=0 no lumping, =1 only keep diagonal entry, =2 accumulate full row to diagonal)"),
    :params => (nothing, "array of parameters that should be made available in qpinfo argument of kernel function"),
    :entry_tolerance => (0, "threshold to add entry to sparse matrix"),
    :use_sparsity_pattern => ("auto", "read sparsity pattern of jacobian of kernel to find out which components couple"),
    :parallel_groups => (false, "assemble operator in parallel using CellAssemblyGroups (assembles separated matrices that are added together sequantially)"),
    :parallel => (false, "assemble operator in parallel using colors/partitions information (assembles into full matrix directly)"),
    :time_dependent => (false, "operator is time-dependent ?"),
    :store => (false, "store matrix separately (and copy from there when reassembly is triggered)"),
    :quadorder => ("auto", "quadrature order"),
    :bonus_quadorder => (0, "additional quadrature order added to quadorder"),
    :verbosity => (0, "verbosity level"),
    :regions => ([], "subset of regions where operator should be assembly only"),
)

getFEStest(FVB::FEVectorBlock) = FVB.FES
getFESansatz(FVB::FEVectorBlock) = FVB.FES
getFEStest(FMB::FEMatrixBlock) = FMB.FES
getFESansatz(FMB::FEMatrixBlock) = FMB.FESY

# informs solver when operator needs reassembly
function depends_nonlinearly_on(O::BilinearOperator)
    return unique(O.u_args)
end

# informs solver in which blocks the operator assembles to
function dependencies_when_linearized(O::BilinearOperator)
    return [unique(O.u_ansatz), unique(O.u_test)]
end

# informs solver when operator needs reassembly in a time dependent setting
function is_timedependent(O::BilinearOperator)
    return O.parameters[:time_dependent]
end

function Base.show(io::IO, O::BilinearOperator)
    dependencies = dependencies_when_linearized(O)
    print(io, "$(O.parameters[:name])($([ansatz_function(dependencies[1][j]) for j in 1:length(dependencies[1])]), $([test_function(dependencies[2][j]) for j in 1:length(dependencies[2])]); entities = $(O.parameters[:entities]))")
    return nothing
end


"""
$(TYPEDSIGNATURES)

Constructs a bilinear operator from a user-supplied matrix `A`, which can be a sparse matrix or a block-structured FEMatrix. The arguments `u_test` and `u_ansatz` specify the test and ansatz (trial) unknowns or indices, determining where the (blocks of the) matrix are inserted in the global system. If `u_ansatz` is not provided, it defaults to `u_test`.

# Arguments
- `A::AbstractMatrix`: The matrix representing the bilinear form, e.g., a sparse matrix.
- `u_test::Vector{<:Union{Unknown, Int}}`: Identifiers or indices for the test functions.
- `u_ansatz::Vector{<:Union{Unknown, Int}}` (optional): Identifiers or indices for the ansatz (trial) functions. Defaults to `u_test`.
- `u_args::Vector{<:Union{Unknown, Int}}` (optional): Identifiers or indices for unknowns on which the matrix depends in a nonlinear way (this tells the solver which solution blocks trigger reassembly).

# Returns
A `BilinearOperatorFromMatrix` object that can be used in the assembly process.


"""
function BilinearOperator(A::AbstractMatrix, u_test::Vector{<:Union{Unknown, Int}}, u_ansatz = u_test, u_args = []; kwargs...)
    parameters = Dict{Symbol, Any}(k => v[1] for (k, v) in default_blfop_kwargs())
    _update_params!(parameters, kwargs)
    return BilinearOperatorFromMatrix{typeof(u_test[1]), typeof(A)}(u_test, u_ansatz, u_args, A, parameters)
end


function BilinearOperator(kernel, u_test, ops_test, u_ansatz = u_test, ops_ansatz = ops_test; Tv = Float64, kwargs...)
    parameters = Dict{Symbol, Any}(k => v[1] for (k, v) in default_blfop_kwargs())
    _update_params!(parameters, kwargs)
    @assert length(u_ansatz) == length(ops_ansatz)
    @assert length(u_test) == length(ops_test)
    if parameters[:store]
        if parameters[:parallel]
            storage = MTExtendableSparseMatrixCSC{Float64, Int}(0, 0, 1)
        else
            storage = ExtendableSparseMatrix{Float64, Int}(0, 0)
        end
    else
        storage = nothing
    end
    return BilinearOperator{Tv, typeof(u_test[1]), typeof(kernel), typeof(storage)}(
        u_test,
        ops_test,
        u_ansatz,
        ops_ansatz,
        [],
        [],
        kernel,
        [[zeros(Tv, 0, 0, 0)]],
        [[zeros(Tv, 0, 0, 0)]],
        [[zeros(Tv, 0, 0, 0)]],
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        storage,
        parameters,
    )
end

function BilinearOperator(kernel, u_test, ops_test, u_ansatz, ops_ansatz, u_args, ops_args; Tv = Float64, kwargs...)
    parameters = Dict{Symbol, Any}(k => v[1] for (k, v) in default_blfop_kwargs())
    _update_params!(parameters, kwargs)
    @assert length(u_args) == length(ops_args)
    @assert length(u_ansatz) == length(ops_ansatz)
    @assert length(u_test) == length(ops_test)
    if parameters[:store]
        storage = ExtendableSparseMatrix{Float64, Int}(0, 0)
    else
        storage = nothing
    end
    return BilinearOperator{Tv, typeof(u_test[1]), typeof(kernel), typeof(storage)}(
        u_test,
        ops_test,
        u_ansatz,
        ops_ansatz,
        u_args,
        ops_args,
        kernel,
        [[zeros(Tv, 0, 0, 0)]],
        [[zeros(Tv, 0, 0, 0)]],
        [[zeros(Tv, 0, 0, 0)]],
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        storage,
        parameters,
    )
end

function BilinearOperator(kernel, oa_test::Array{<:Tuple{Union{Unknown, Int}, DataType}, 1}, oa_ansatz::Array{<:Tuple{Union{Unknown, Int}, DataType}, 1} = oa_test; kwargs...)
    u_test = [oa[1] for oa in oa_test]
    u_ansatz = [oa[1] for oa in oa_ansatz]
    ops_test = [oa[2] for oa in oa_test]
    ops_ansatz = [oa[2] for oa in oa_ansatz]
    return BilinearOperator(kernel, u_test, ops_test, u_ansatz, ops_ansatz; kwargs...)
end


"""
````
function BilinearOperator(
	[kernel!::Function],
	oa_test::Array{<:Tuple{Union{Unknown,Int}, DataType},1},
	oa_ansatz::Array{<:Tuple{Union{Unknown,Int}, DataType},1} = oa_test;
	kwargs...)
````

Constructs a bilinear form by evaluating the vector product of operator evaluations of the test and ansatz functions. Operator evaluations are specified as tuples pairing an unknown identifier (or integer) with a function operator (such as those generated by `grad(u)`, `id(u)`, etc).

If a kernel function is provided as the first argument, it customizes the ansatz function evaluation. The kernel must have the signature:

    kernel!(result, eval_ansatz, qpinfo)

where `result` is the output vector, `eval_ansatz` is the ansatz function evaluation at the quadrature point, and `qpinfo` provides quadrature point information.

If no kernel is provided, the standard kernel is used that produces the dot product.

# Arguments
- `oa_test`: Array of tuples `(unknown, operator)` for test functions.
- `oa_ansatz`: Array of tuples `(unknown, operator)` for ansatz (trial) functions. Defaults to `oa_test`.

# Keyword Arguments
$(_myprint(default_blfop_kwargs()))

# Example
```julia
BilinearOperator([grad(1)], [grad(1)]; entities=ON_CELLS)
```

"""
function BilinearOperator(oa_test::Array{<:Tuple{Union{Unknown, Int}, DataType}, 1}, oa_ansatz::Array{<:Tuple{Union{Unknown, Int}, DataType}, 1} = oa_test; kwargs...)
    u_test = [oa[1] for oa in oa_test]
    u_ansatz = [oa[1] for oa in oa_ansatz]
    ops_test = [oa[2] for oa in oa_test]
    ops_ansatz = [oa[2] for oa in oa_ansatz]
    return BilinearOperator(ExtendableFEMBase.standard_kernel, u_test, ops_test, u_ansatz, ops_ansatz; kwargs...)
end


"""
$(TYPEDSIGNATURES)

Constructs a nonlinear bilinear form by evaluating a user-supplied kernel function that depends on both the operator evaluations of the ansatz (trial) functions and additional argument functions (e.g., the current solution). The result of the kernel is then contracted with the operator evaluations of the test functions. This is typically used for linearizations of nonlinear operators.

The kernel function must have the following signature:

    kernel!(result, eval_ansatz, eval_args, qpinfo)

where
- `result`: Output vector for the kernel evaluation.
- `eval_ansatz`: Operator evaluation(s) of the ansatz functions at the quadrature point.
- `eval_args`: Operator evaluation(s) of the argument functions (e.g., current solution) at the quadrature point.
- `qpinfo`: Quadrature point information structure.

Operator evaluations are specified as tuples pairing an unknown identifier (or integer) with a function operator (such as those generated by `grad(u)`, `id(u)`, etc).

# Arguments
- `kernel`: User-supplied kernel function as described above.
- `oa_test`: Array of tuples `(unknown, operator)` for test functions.
- `oa_ansatz`: Array of tuples `(unknown, operator)` for ansatz (trial) functions.
- `oa_args`: Array of tuples `(unknown, operator)` for argument functions (e.g., current solution).
- `kwargs...`: Additional keyword arguments for operator assembly (see below).

# Keyword Arguments
$(_myprint(default_blfop_kwargs()))

"""
function BilinearOperator(kernel, oa_test::Array{<:Tuple{Union{Unknown, Int}, DataType}, 1}, oa_ansatz::Array{<:Tuple{Union{Unknown, Int}, DataType}, 1}, oa_args::Array{<:Tuple{Union{Unknown, Int}, DataType}, 1}; kwargs...)
    u_test = [oa[1] for oa in oa_test]
    u_ansatz = [oa[1] for oa in oa_ansatz]
    u_args = [oa[1] for oa in oa_args]
    ops_test = [oa[2] for oa in oa_test]
    ops_ansatz = [oa[2] for oa in oa_ansatz]
    ops_args = [oa[2] for oa in oa_args]
    return BilinearOperator(kernel, u_test, ops_test, u_ansatz, ops_ansatz, u_args, ops_args; kwargs...)
end

function build_assembler!(A::AbstractMatrix, O::BilinearOperator{Tv}, FE_test, FE_ansatz, FE_args::Array{<:FEVectorBlock, 1}; time = 0.0, kwargs...) where {Tv}
    ## check if FES is the same as last time
    FES_test = [getFEStest(FE_test[j]) for j in 1:length(FE_test)]
    FES_ansatz = [getFESansatz(FE_ansatz[j]) for j in 1:length(FE_ansatz)]
    FES_args = [FE_args[j].FES for j in 1:length(FE_args)]
    return if (O.FES_test != FES_test) || (O.FES_args != FES_args)

        if O.parameters[:verbosity] > 0
            @info "$(O.parameters[:name]) : building assembler"
        end

        ## determine grid
        xgrid = determine_assembly_grid(FES_test, FES_ansatz)

        ## prepare assembly
        AT = O.parameters[:entities]
        Ti = typeof(xgrid).parameters[2]
        if xgrid == FES_test[1].dofgrid
            gridAT = ExtendableFEMBase.EffAT4AssemblyType(get_AT(FES_test[1]), AT)
        else
            gridAT = AT
        end

        itemgeometries = xgrid[GridComponentGeometries4AssemblyType(gridAT)]
        itemvolumes = xgrid[GridComponentVolumes4AssemblyType(gridAT)]
        itemregions = xgrid[GridComponentRegions4AssemblyType(gridAT)]
        if num_pcolors(xgrid) > 1 && gridAT == ON_CELLS
            maxnpartitions = maximum(num_partitions_per_color(xgrid))
            pc = xgrid[PartitionCells]
            itemassemblygroups = [pc[j]:(pc[j + 1] - 1) for j in 1:num_partitions(xgrid)]
            # assuming here that all cells of one partition have the same geometry
        else
            itemassemblygroups = xgrid[GridComponentAssemblyGroups4AssemblyType(gridAT)]
            itemassemblygroups = [view(itemassemblygroups, :, j) for j in 1:num_sources(itemassemblygroups)]
        end
        has_normals = true
        if AT <: ON_FACES
            itemnormals = xgrid[FaceNormals]
        elseif AT <: ON_BFACES
            itemnormals = xgrid[FaceNormals][:, xgrid[BFaceFaces]]
        else
            has_normals = false
        end
        FETypes_test = [eltype(F) for F in FES_test]
        FETypes_ansatz = [eltype(F) for F in FES_ansatz]
        FETypes_args = [eltype(F) for F in FES_args]
        EGs = [itemgeometries[itemassemblygroups[j][1]] for j in 1:length(itemassemblygroups)]

        ## prepare assembly
        nargs = length(FES_args)
        ntest = length(FES_test)
        nansatz = length(FES_ansatz)
        O.QF = []
        O.BE_test = Array{Array{<:FEEvaluator, 1}, 1}([])
        O.BE_ansatz = Array{Array{<:FEEvaluator, 1}, 1}([])
        O.BE_args = Array{Array{<:FEEvaluator, 1}, 1}([])
        O.BE_test_vals = Array{Array{Array{Tv, 3}, 1}, 1}([])
        O.BE_ansatz_vals = Array{Array{Array{Tv, 3}, 1}, 1}([])
        O.BE_args_vals = Array{Array{Array{Tv, 3}, 1}, 1}([])
        O.QP_infos = Array{QPInfos, 1}([])
        O.L2G = []
        for EG in EGs
            ## quadrature formula for EG
            polyorder_ansatz = maximum([get_polynomialorder(FETypes_ansatz[j], EG) - ExtendableFEMBase.NeededDerivative4Operator(O.ops_ansatz[j]) for j in 1:nansatz])
            polyorder_test = maximum([get_polynomialorder(FETypes_test[j], EG) - ExtendableFEMBase.NeededDerivative4Operator(O.ops_test[j]) for j in 1:ntest])
            if O.parameters[:quadorder] == "auto"
                quadorder = polyorder_ansatz + polyorder_test + O.parameters[:bonus_quadorder]
            else
                quadorder = O.parameters[:quadorder] + O.parameters[:bonus_quadorder]
            end
            if O.parameters[:verbosity] > 1
                @info "...... integrating on $EG with quadrature order $quadorder"
            end
            push!(O.QF, QuadratureRule{Tv, EG}(quadorder))

            ## L2G map for EG
            push!(O.L2G, L2GTransformer(EG, xgrid, gridAT))

            ## FE basis evaluator for EG
            push!(O.BE_test, [FEEvaluator(FES_test[j], O.ops_test[j], O.QF[end], xgrid; AT = AT, L2G = O.L2G[end]) for j in 1:ntest])
            push!(O.BE_ansatz, [FEEvaluator(FES_ansatz[j], O.ops_ansatz[j], O.QF[end], xgrid; AT = AT, L2G = O.L2G[end]) for j in 1:nansatz])
            push!(O.BE_args, [FEEvaluator(FES_args[j], O.ops_args[j], O.QF[end], xgrid; AT = AT, L2G = O.L2G[end]) for j in 1:nargs])
            push!(O.BE_test_vals, [BE.cvals for BE in O.BE_test[end]])
            push!(O.BE_ansatz_vals, [BE.cvals for BE in O.BE_ansatz[end]])
            push!(O.BE_args_vals, [BE.cvals for BE in O.BE_args[end]])

            ## parameter structure
            push!(O.QP_infos, QPInfos(xgrid; time = time, params = O.parameters[:params]))
        end

        ## prepare regions
        regions = O.parameters[:regions]
        visit_region = zeros(Bool, maximum(itemregions))
        if length(regions) > 0
            visit_region[regions] .= true
        else
            visit_region .= true
        end

        ## prepare operator infos
        op_lengths_test = [size(O.BE_test[1][j].cvals, 1) for j in 1:ntest]
        op_lengths_ansatz = [size(O.BE_ansatz[1][j].cvals, 1) for j in 1:nansatz]
        op_lengths_args = [size(O.BE_args[1][j].cvals, 1) for j in 1:nargs]

        op_offsets_test = [0]
        op_offsets_ansatz = [0]
        op_offsets_args = [0]
        append!(op_offsets_test, cumsum(op_lengths_test))
        append!(op_offsets_ansatz, cumsum(op_lengths_ansatz))
        append!(op_offsets_args, cumsum(op_lengths_args))
        offsets_test = [FE_test[j].offset for j in 1:length(FES_test)]
        offsets_ansatz = [FE_ansatz[j].offsetY for j in 1:length(FES_ansatz)]

        ## prepare sparsity pattern
        use_sparsity_pattern = O.parameters[:use_sparsity_pattern]
        if use_sparsity_pattern == "auto"
            use_sparsity_pattern = false
        end
        coupling_matrix::Matrix{Bool} = ones(Bool, nansatz, ntest)
        if use_sparsity_pattern
            kernel_params = (result, input) -> (O.kernel(result, input, O.QP_infos[1]))
            detector = TracerSparsityDetector()
            sparsity_pattern = jacobian_sparsity(kernel_params, zeros(Tv, op_offsets_test[end]), zeros(Tv, op_offsets_ansatz[end]), detector)

            ## find out which test and ansatz functions couple
            for id in 1:nansatz
                for idt in 1:ntest
                    couple = false
                    for j in 1:op_lengths_ansatz[id]
                        for k in 1:op_lengths_test[idt]
                            if sparsity_pattern[k + op_offsets_test[idt], j + op_offsets_ansatz[id]] > 0
                                couple = true
                            end
                        end
                    end
                    coupling_matrix[id, idt] = couple
                end
            end
        end
        couples_with::Vector{Vector{Int}} = [findall(==(true), view(coupling_matrix, j, :)) for j in 1:nansatz]

        ## prepare parallel assembly_allocations
        if O.parameters[:parallel_groups]
            Aj = Array{typeof(A), 1}(undef, length(EGs))
            for j in 1:length(EGs)
                Aj[j] = copy(A)
            end
        end

        FEATs_test = [ExtendableFEMBase.EffAT4AssemblyType(get_AT(FES_test[j]), AT) for j in 1:ntest]
        FEATs_ansatz = [ExtendableFEMBase.EffAT4AssemblyType(get_AT(FES_ansatz[j]), AT) for j in 1:nansatz]
        FEATs_args = [ExtendableFEMBase.EffAT4AssemblyType(get_AT(FES_args[j]), AT) for j in 1:nargs]
        itemdofs_test::Array{Union{Adjacency{Ti}, SerialVariableTargetAdjacency{Ti}}, 1} = [get_dofmap(FES_test[j], xgrid, FEATs_test[j]) for j in 1:ntest]
        itemdofs_ansatz::Array{Union{Adjacency{Ti}, SerialVariableTargetAdjacency{Ti}}, 1} = [get_dofmap(FES_ansatz[j], xgrid, FEATs_ansatz[j]) for j in 1:nansatz]
        itemdofs_args::Array{Union{Adjacency{Ti}, SerialVariableTargetAdjacency{Ti}}, 1} = [get_dofmap(FES_args[j], xgrid, FEATs_args[j]) for j in 1:nargs]
        factor = O.parameters[:factor]
        transposed_copy = O.parameters[:transposed_copy]
        entry_tol = O.parameters[:entry_tolerance]
        lump = O.parameters[:lump]

        ## Assembly loop for fixed geometry
        function assembly_loop(
                A::AbstractSparseArray{T},
                sol::Array{<:FEVectorBlock, 1},
                items,
                EG::ElementGeometries,
                QF::QuadratureRule,
                BE_test::Array{<:FEEvaluator, 1},
                BE_ansatz::Array{<:FEEvaluator, 1},
                BE_args::Array{<:FEEvaluator, 1},
                BE_test_vals::Array{Array{Tv, 3}, 1},
                BE_ansatz_vals::Array{Array{Tv, 3}, 1},
                BE_args_vals::Array{Array{Tv, 3}, 1},
                L2G::L2GTransformer,
                QPinfos::QPInfos,
                part = 1,
            ) where {T}

            input_ansatz = zeros(T, op_offsets_ansatz[end])
            input_args = zeros(T, op_offsets_args[end])
            result_kernel = zeros(T, op_offsets_test[end])

            ndofs_test::Array{Int, 1} = [size(BE.cvals, 2) for BE in BE_test]
            ndofs_ansatz::Array{Int, 1} = [size(BE.cvals, 2) for BE in BE_ansatz]
            ndofs_args::Array{Int, 1} = [size(BE.cvals, 2) for BE in BE_args]

            Aloc = Matrix{Matrix{T}}(undef, ntest, nansatz)
            for j in 1:ntest, k in 1:nansatz
                Aloc[j, k] = zeros(T, ndofs_test[j], ndofs_ansatz[k])
            end
            weights, xref = QF.w, QF.xref
            nweights = length(weights)

            for item::Int in items
                if itemregions[item] > 0
                    if !(visit_region[itemregions[item]]) || AT == ON_IFACES
                        continue
                    end
                end
                QPinfos.region = itemregions[item]
                QPinfos.item = item
                if has_normals
                    QPinfos.normal .= view(itemnormals, :, item)
                end
                QPinfos.volume = itemvolumes[item]

                ## update FE basis evaluators
                for j in 1:ntest
                    BE_test[j].citem[] = item
                    update_basis!(BE_test[j])
                end
                for j in 1:nansatz
                    BE_ansatz[j].citem[] = item
                    update_basis!(BE_ansatz[j])
                end
                for j in 1:nargs
                    BE_args[j].citem[] = item
                    update_basis!(BE_args[j])
                end
                update_trafo!(L2G, item)

                ## evaluate arguments
                for qp in 1:nweights
                    fill!(input_args, 0)
                    for id in 1:nargs
                        for j in 1:ndofs_args[id]
                            dof_j = itemdofs_args[id][j, item]
                            for d in 1:op_lengths_args[id]
                                input_args[d + op_offsets_args[id]] += sol[id][dof_j] * BE_args_vals[id][d, j, qp]
                            end
                        end
                    end

                    ## get global x for quadrature point
                    eval_trafo!(QPinfos.x, L2G, xref[qp])

                    # update matrix
                    for id in 1:nansatz
                        for j in 1:ndofs_ansatz[id]
                            # evaluat kernel for ansatz basis function
                            fill!(input_ansatz, 0)
                            for d in 1:op_lengths_ansatz[id]
                                input_ansatz[d + op_offsets_ansatz[id]] += BE_ansatz_vals[id][d, j, qp]
                            end

                            # evaluate kernel
                            O.kernel(result_kernel, input_ansatz, input_args, QPinfos)
                            result_kernel .*= factor * weights[qp]

                            # multiply test function operator evaluation
                            if lump == 1
                                for d in 1:op_lengths_test[id]
                                    Aloc[id, id][j, j] += result_kernel[d + op_offsets_test[id]] * BE_test_vals[id][d, j, qp]
                                end
                            elseif lump == 2
                                for k in 1:ndofs_test[id]
                                    for d in 1:op_lengths_test[id]
                                        Aloc[id, id][j, j] += result_kernel[d + op_offsets_test[id]] * BE_test_vals[id][d, k, qp]
                                    end
                                end
                            else
                                for idt in couples_with[id]
                                    for k in 1:ndofs_test[idt]
                                        for d in 1:op_lengths_test[idt]
                                            Aloc[idt, id][k, j] += result_kernel[d + op_offsets_test[idt]] * BE_test_vals[idt][d, k, qp]
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                ## add local matrices to global matrix
                for id in 1:nansatz, idt in 1:ntest
                    Aloc[idt, id] .*= itemvolumes[item]
                    for j in 1:ndofs_test[idt]
                        dof_j = itemdofs_test[idt][j, item] + offsets_test[idt]
                        for k in 1:ndofs_ansatz[id]
                            dof_k = itemdofs_ansatz[id][k, item] + offsets_ansatz[id]
                            if abs(Aloc[idt, id][j, k]) > entry_tol
                                rawupdateindex!(A, +, Aloc[idt, id][j, k], dof_j, dof_k, part)
                            end
                        end
                    end
                end
                if transposed_copy != 0
                    for id in 1:nansatz, idt in 1:ntest
                        Aloc[idt, id] .*= transposed_copy
                        for j in 1:ndofs_test[idt]
                            dof_j = itemdofs_test[idt][j, item] + offsets_test[idt]
                            for k in 1:ndofs_ansatz[id]
                                dof_k = itemdofs_ansatz[id][k, item] + offsets_ansatz[id]
                                if abs(Aloc[idt, id][j, k]) > entry_tol
                                    rawupdateindex!(A, +, Aloc[idt, id][j, k], dof_k, dof_j, part)
                                end
                            end
                        end
                    end
                end

                for id in 1:nansatz, idt in 1:ntest
                    fill!(Aloc[idt, id], 0)
                end
            end
            return
        end
        O.FES_test = FES_test
        O.FES_ansatz = FES_ansatz
        O.FES_args = FES_args

        function assembler(A, b, sol; kwargs...)
            time = @elapsed begin
                if O.parameters[:parallel]
                    pcp = xgrid[PColorPartitions]
                    ncolors = length(pcp) - 1
                    if O.parameters[:verbosity] > 0
                        @info "$(O.parameters[:name]) : assembling in parallel with $ncolors colors, $(length(EGs)) partitions and $(Threads.nthreads()) threads"
                    end
                    for color in 1:ncolors
                        Threads.@threads for part in pcp[color]:(pcp[color + 1] - 1)
                            assembly_loop(
                                A,
                                sol,
                                itemassemblygroups[part],
                                EGs[part],
                                O.QF[part],
                                O.BE_test[part],
                                O.BE_ansatz[part],
                                O.BE_args[part],
                                O.BE_test_vals[part],
                                O.BE_ansatz_vals[part],
                                O.BE_args_vals[part],
                                O.L2G[part],
                                O.QP_infos[part],
                                part;
                                kwargs...,
                            )
                        end
                    end
                elseif O.parameters[:parallel_groups]
                    Threads.@threads for j in 1:length(EGs)
                        fill!(Aj[j].cscmatrix.nzval, 0)
                        assembly_loop(Aj[j], sol, itemassemblygroups[j], EGs[j], O.QF[j], O.BE_test[j], O.BE_ansatz[j], O.BE_args[j], O.BE_test_vals[j], O.BE_ansatz_vals[j], O.BE_args_vals[j], O.L2G[j], O.QP_infos[j]; kwargs...)
                        flush!(Aj[j])
                    end
                    for j in 1:length(EGs)
                        add!(A, Aj[j])
                    end
                else
                    for j in 1:length(EGs)
                        assembly_loop(A, sol, itemassemblygroups[j], EGs[j], O.QF[j], O.BE_test[j], O.BE_ansatz[j], O.BE_args[j], O.BE_test_vals[j], O.BE_ansatz_vals[j], O.BE_args_vals[j], O.L2G[j], O.QP_infos[j]; kwargs...)
                    end
                end
            end
            flush!(A)
            return if O.parameters[:verbosity] > 0
                @info "$(O.parameters[:name]) : assembly took $time s"
            end
        end
        O.assembler = assembler
    else
        ## update the time
        for j in 1:length(O.QP_infos)
            O.QP_infos[j].time = time
        end
    end
end

function build_assembler!(A, O::BilinearOperator{Tv}, FE_test, FE_ansatz; time = 0.0, kwargs...) where {Tv}
    ## check if FES is the same as last time
    FES_test = [getFEStest(FE_test[j]) for j in 1:length(FE_test)]
    FES_ansatz = [getFESansatz(FE_ansatz[j]) for j in 1:length(FE_ansatz)]

    return if (O.FES_test != FES_test) || (O.FES_ansatz != FES_ansatz)

        ntest = length(FES_test)
        nansatz = length(FES_ansatz)
        if O.parameters[:verbosity] > 0
            @info "$(O.parameters[:name]) : building assembler"
            if O.parameters[:verbosity] > 1
                @info "......   TEST : $([(get_FEType(FES_test[j]), O.ops_test[j]) for j in 1:ntest])"
                @info "...... ANSATZ : $([(get_FEType(FES_ansatz[j]), O.ops_ansatz[j]) for j in 1:nansatz])"
            end
        end

        ## determine grid
        xgrid = determine_assembly_grid(FES_test, FES_ansatz)

        ## prepare assembly
        AT = O.parameters[:entities]
        Ti = typeof(xgrid).parameters[2]
        if xgrid == FES_test[1].dofgrid
            gridAT = ExtendableFEMBase.EffAT4AssemblyType(get_AT(FES_test[1]), AT)
        else
            gridAT = AT
        end
        if O.parameters[:verbosity] > 1
            @info "......     AT : $(AT)"
            @info "...... gridAT : $(gridAT)"
        end


        itemgeometries = xgrid[GridComponentGeometries4AssemblyType(gridAT)]
        itemvolumes = xgrid[GridComponentVolumes4AssemblyType(gridAT)]
        itemregions = xgrid[GridComponentRegions4AssemblyType(gridAT)]
        if num_pcolors(xgrid) > 1 && gridAT == ON_CELLS
            maxnpartitions = maximum(num_partitions_per_color(xgrid))
            pc = xgrid[PartitionCells]
            itemassemblygroups = [pc[j]:(pc[j + 1] - 1) for j in 1:num_partitions(xgrid)]
            # assuming here that all cells of one partition have the same geometry
        else
            itemassemblygroups = xgrid[GridComponentAssemblyGroups4AssemblyType(gridAT)]
            itemassemblygroups = [view(itemassemblygroups, :, j) for j in 1:num_sources(itemassemblygroups)]
        end
        EGs = [itemgeometries[itemassemblygroups[j][1]] for j in 1:length(itemassemblygroups)]

        has_normals = true
        if gridAT <: ON_FACES
            itemnormals = xgrid[FaceNormals]
        elseif gridAT <: ON_BFACES
            itemnormals = xgrid[FaceNormals][:, xgrid[BFaceFaces]]
        else
            has_normals = false
        end
        FETypes_test = [eltype(F) for F in FES_test]
        FETypes_ansatz = [eltype(F) for F in FES_ansatz]

        ## prepare assembly
        O.QF = []
        O.BE_test = Array{Array{<:FEEvaluator, 1}, 1}([])
        O.BE_ansatz = Array{Array{<:FEEvaluator, 1}, 1}([])
        O.BE_test_vals = Array{Array{Array{Tv, 3}, 1}, 1}([])
        O.BE_ansatz_vals = Array{Array{Array{Tv, 3}, 1}, 1}([])
        O.QP_infos = Array{QPInfos, 1}([])
        O.L2G = []
        for EG in EGs
            ## quadrature formula for EG
            polyorder_ansatz = maximum([get_polynomialorder(FETypes_ansatz[j], EG) - ExtendableFEMBase.NeededDerivative4Operator(O.ops_ansatz[j]) for j in 1:nansatz])
            polyorder_test = maximum([get_polynomialorder(FETypes_test[j], EG) - ExtendableFEMBase.NeededDerivative4Operator(O.ops_test[j]) for j in 1:ntest])
            if O.parameters[:quadorder] == "auto"
                quadorder = polyorder_ansatz + polyorder_test + O.parameters[:bonus_quadorder]
            else
                quadorder = O.parameters[:quadorder] + O.parameters[:bonus_quadorder]
            end
            if O.parameters[:verbosity] > 1
                @info "...... integrating on $EG with quadrature order $quadorder"
            end
            push!(O.QF, QuadratureRule{Tv, EG}(quadorder))

            ## L2G map for EG
            push!(O.L2G, L2GTransformer(EG, xgrid, gridAT))

            ## FE basis evaluator for EG
            push!(O.BE_test, [FEEvaluator(FES_test[j], O.ops_test[j], O.QF[end], xgrid; AT = AT, L2G = O.L2G[end]) for j in 1:ntest])
            push!(O.BE_ansatz, [FEEvaluator(FES_ansatz[j], O.ops_ansatz[j], O.QF[end], xgrid; AT = AT, L2G = O.L2G[end]) for j in 1:nansatz])
            push!(O.BE_test_vals, [BE.cvals for BE in O.BE_test[end]])
            push!(O.BE_ansatz_vals, [BE.cvals for BE in O.BE_ansatz[end]])

            ## parameter structure
            push!(O.QP_infos, QPInfos(xgrid; time = time, x = ones(Tv, size(xgrid[Coordinates], 1)), params = O.parameters[:params]))
        end

        ## prepare regions
        regions = O.parameters[:regions]
        visit_region = zeros(Bool, maximum(itemregions))
        if length(regions) > 0
            visit_region[regions] .= true
        else
            visit_region .= true
        end

        ## prepare operator infos
        op_lengths_test = [size(O.BE_test[1][j].cvals, 1) for j in 1:ntest]
        op_lengths_ansatz = [size(O.BE_ansatz[1][j].cvals, 1) for j in 1:nansatz]

        op_offsets_test = [0]
        op_offsets_ansatz = [0]
        append!(op_offsets_test, cumsum(op_lengths_test))
        append!(op_offsets_ansatz, cumsum(op_lengths_ansatz))
        offsets_test = [FE_test[j].offset for j in 1:length(FES_test)]
        offsets_ansatz = [FE_ansatz[j].offsetY for j in 1:length(FES_ansatz)]

        ## prepare sparsity pattern
        use_sparsity_pattern = O.parameters[:use_sparsity_pattern]
        if use_sparsity_pattern == "auto"
            use_sparsity_pattern = ntest > 1
        end
        coupling_matrix::Matrix{Bool} = ones(Bool, nansatz, ntest)
        if use_sparsity_pattern
            kernel_params = (result, input) -> (O.kernel(result, input, O.QP_infos[1]))
            detector = TracerSparsityDetector()
            sparsity_pattern = jacobian_sparsity(kernel_params, zeros(Tv, op_offsets_test[end]), zeros(Tv, op_offsets_ansatz[end]), detector)

            ## find out which test and ansatz functions couple
            for id in 1:nansatz
                for idt in 1:ntest
                    couple = false
                    for j in 1:op_lengths_ansatz[id]
                        for k in 1:op_lengths_test[idt]
                            if sparsity_pattern[k + op_offsets_test[idt], j + op_offsets_ansatz[id]] > 0
                                couple = true
                            end
                        end
                    end
                    coupling_matrix[id, idt] = couple
                end
            end
        end
        couples_with::Vector{Vector{Int}} = [findall(==(true), view(coupling_matrix, j, :)) for j in 1:nansatz]

        ## prepare parallel assembly
        if O.parameters[:parallel_groups]
            Aj = Array{typeof(A), 1}(undef, length(EGs))
            for j in 1:length(EGs)
                Aj[j] = deepcopy(A)
            end
        end

        FEATs_test = [ExtendableFEMBase.EffAT4AssemblyType(get_AT(FES_test[j]), AT) for j in 1:ntest]
        FEATs_ansatz = [ExtendableFEMBase.EffAT4AssemblyType(get_AT(FES_ansatz[j]), AT) for j in 1:nansatz]

        itemdofs_test::Array{Union{Adjacency{Ti}, SerialVariableTargetAdjacency{Ti}}, 1} = [get_dofmap(FES_test[j], xgrid, FEATs_test[j]) for j in 1:ntest]
        itemdofs_ansatz::Array{Union{Adjacency{Ti}, SerialVariableTargetAdjacency{Ti}}, 1} = [get_dofmap(FES_ansatz[j], xgrid, FEATs_ansatz[j]) for j in 1:nansatz]
        factor = O.parameters[:factor]
        transposed_copy = O.parameters[:transposed_copy]
        entry_tol = O.parameters[:entry_tolerance]
        lump = O.parameters[:lump]

        ## Assembly loop for fixed geometry
        function assembly_loop(
                A::AbstractSparseArray{T},
                items,
                EG::ElementGeometries,
                QF::QuadratureRule,
                BE_test::Array{<:FEEvaluator, 1},
                BE_ansatz::Array{<:FEEvaluator, 1},
                BE_test_vals::Array{Array{Tv, 3}, 1},
                BE_ansatz_vals::Array{Array{Tv, 3}, 1},
                L2G::L2GTransformer,
                QPinfos::QPInfos,
                part = 1,
            ) where {T}

            input_ansatz = zeros(T, op_offsets_ansatz[end])
            result_kernel = zeros(T, op_offsets_test[end])

            ndofs_test::Array{Int, 1} = [size(BE.cvals, 2) for BE in BE_test]
            ndofs_ansatz::Array{Int, 1} = [size(BE.cvals, 2) for BE in BE_ansatz]

            Aloc = Matrix{Matrix{T}}(undef, ntest, nansatz)
            for j in 1:ntest, k in 1:nansatz
                Aloc[j, k] = zeros(T, ndofs_test[j], ndofs_ansatz[k])
            end
            weights, xref = QF.w, QF.xref
            nweights = length(weights)

            for item::Int in items
                if itemregions[item] > 0
                    if !(visit_region[itemregions[item]]) || AT == ON_IFACES
                        continue
                    end
                else
                    if length(regions) > 0
                        continue
                    end
                end
                QPinfos.region = itemregions[item]
                QPinfos.item = item
                if has_normals
                    QPinfos.normal .= view(itemnormals, :, item)
                end
                QPinfos.volume = itemvolumes[item]

                ## update FE basis evaluators
                for j in 1:ntest
                    BE_test[j].citem[] = item
                    update_basis!(BE_test[j])
                end
                for j in 1:nansatz
                    BE_ansatz[j].citem[] = item
                    update_basis!(BE_ansatz[j])
                end
                update_trafo!(L2G, item)

                ## evaluate arguments
                for qp in 1:nweights

                    ## get global x for quadrature point
                    eval_trafo!(QPinfos.x, L2G, xref[qp])

                    # update matrix
                    for id in 1:nansatz
                        for j in 1:ndofs_ansatz[id]
                            # evaluat kernel for ansatz basis function
                            fill!(input_ansatz, 0)
                            for d in 1:op_lengths_ansatz[id]
                                input_ansatz[d + op_offsets_ansatz[id]] = BE_ansatz_vals[id][d, j, qp]
                            end

                            # evaluate kernel
                            O.kernel(result_kernel, input_ansatz, QPinfos)
                            result_kernel .*= factor * weights[qp]

                            # multiply test function operator evaluation
                            if lump == 1
                                for d in 1:op_lengths_test[id]
                                    Aloc[id, id][j, j] += result_kernel[d + op_offsets_test[id]] * BE_test_vals[id][d, j, qp]
                                end
                            elseif lump == 2
                                for k in 1:ndofs_test[id]
                                    for d in 1:op_lengths_test[id]
                                        Aloc[id, id][j, j] += result_kernel[d + op_offsets_test[id]] * BE_test_vals[id][d, k, qp]
                                    end
                                end
                            else
                                for idt in couples_with[id]
                                    for k in 1:ndofs_test[idt]
                                        for d in 1:op_lengths_test[idt]
                                            Aloc[idt, id][k, j] += result_kernel[d + op_offsets_test[idt]] * BE_test_vals[idt][d, k, qp]
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                ## add local matrices to global matrix
                for id in 1:nansatz, idt in couples_with[id]
                    Aloc[idt, id] .*= itemvolumes[item]
                    for j in 1:ndofs_test[idt]
                        dof_j = itemdofs_test[idt][j, item] + offsets_test[idt]
                        for k in 1:ndofs_ansatz[id]
                            dof_k = itemdofs_ansatz[id][k, item] + offsets_ansatz[id]
                            if abs(Aloc[idt, id][j, k]) > entry_tol
                                rawupdateindex!(A, +, Aloc[idt, id][j, k], dof_j, dof_k, part)
                            end
                        end
                    end
                end
                if transposed_copy != 0
                    for id in 1:nansatz, idt in couples_with[id]
                        Aloc[idt, id] .*= transposed_copy
                        for j in 1:ndofs_test[idt]
                            dof_j = itemdofs_test[idt][j, item] + offsets_test[idt]
                            for k in 1:ndofs_ansatz[id]
                                dof_k = itemdofs_ansatz[id][k, item] + offsets_ansatz[id]
                                if abs(Aloc[idt, id][j, k]) > entry_tol
                                    rawupdateindex!(A, +, Aloc[idt, id][j, k], dof_k, dof_j, part)
                                end
                            end
                        end
                    end
                end

                for id in 1:nansatz, idt in 1:ntest
                    fill!(Aloc[idt, id], 0)
                end
            end
            return
        end
        O.FES_test = FES_test
        O.FES_ansatz = FES_ansatz

        function assembler(A, b; kwargs...)
            time = @elapsed begin
                if O.parameters[:store] && size(A) == size(O.storage)
                    add!(A, O.storage)
                else
                    if O.parameters[:store]
                        if O.parameters[:parallel]
                            S = MTExtendableSparseMatrixCSC{Float64, Int}(size(A, 1), size(A, 2), length(EGs))
                        else
                            S = ExtendableSparseMatrix{Float64, Int}(size(A, 1), size(A, 2))
                        end
                    else
                        S = A
                    end
                    if O.parameters[:parallel]
                        pcp = xgrid[PColorPartitions]
                        ncolors = length(pcp) - 1
                        if O.parameters[:verbosity] > 0
                            @info "$(O.parameters[:name]) : assembling in parallel with $ncolors colors, $(length(EGs)) partitions and $(Threads.nthreads()) threads"
                        end
                        for color in 1:ncolors
                            Threads.@threads for part in pcp[color]:(pcp[color + 1] - 1)
                                assembly_loop(S, itemassemblygroups[part], EGs[part], O.QF[part], O.BE_test[part], O.BE_ansatz[part], O.BE_test_vals[part], O.BE_ansatz_vals[part], O.L2G[part], O.QP_infos[part], part; kwargs...)
                            end
                        end
                    elseif O.parameters[:parallel_groups]
                        Threads.@threads for j in 1:length(EGs)
                            fill!(Aj[j].cscmatrix.nzval, 0)
                            assembly_loop(Aj[j], itemassemblygroups[j], EGs[j], O.QF[j], O.BE_test[j], O.BE_ansatz[j], O.BE_test_vals[j], O.BE_ansatz_vals[j], O.L2G[j], O.QP_infos[j], j; kwargs...)
                        end
                        for j in 1:length(EGs)
                            add!(S, Aj[j])
                        end
                    else
                        for j in 1:length(EGs)
                            assembly_loop(S, itemassemblygroups[j], EGs[j], O.QF[j], O.BE_test[j], O.BE_ansatz[j], O.BE_test_vals[j], O.BE_ansatz_vals[j], O.L2G[j], O.QP_infos[j]; kwargs...)
                        end
                    end
                    flush!(S)
                    if O.parameters[:store]
                        add!(A, S)
                        O.storage = S
                    end
                end
            end
            return if O.parameters[:verbosity] > 0
                @info "$(O.parameters[:name]) : assembly took $time s"
            end
        end
        O.assembler = assembler
    else
        ## update the time
        for j in 1:length(O.QP_infos)
            O.QP_infos[j].time = time
        end
    end
end

function assemble!(A, b, sol, O::BilinearOperator{Tv, UT}, SC::SolverConfiguration; assemble_matrix = true, assemblr_rhs = true, kwargs...) where {Tv, UT}
    if !assemble_matrix
        return nothing
    end
    if UT <: Integer
        ind_test = O.u_test
        ind_ansatz = O.u_ansatz
        ind_args = O.u_args
    elseif UT <: Unknown
        ind_test = [get_unknown_id(SC, u) for u in O.u_test]
        ind_ansatz = [get_unknown_id(SC, u) for u in O.u_ansatz]
        ind_args = [findfirst(==(u), sol.tags) for u in O.u_args] #[get_unknown_id(SC, u) for u in O.u_args]
    end
    return if length(O.u_args) > 0
        build_assembler!(A.entries, O, [A[j, j] for j in ind_test], [A[j, j] for j in ind_ansatz], [sol[j] for j in ind_args]; kwargs...)
        O.assembler(A.entries, b.entries, [sol[j] for j in ind_args])
    else
        build_assembler!(A.entries, O, [A[j, j] for j in ind_test], [A[j, j] for j in ind_ansatz]; kwargs...)
        O.assembler(A.entries, b.entries)
    end
end

"""
$(TYPEDSIGNATURES)

Assembles the given `BilinearOperator` into the provided `FEMatrix` `A`. This function computes the matrix representation of the bilinear form defined by `O` and inserts the resulting blocks into `A`. If the operator depends on additional argument functions (e.g., for nonlinear problems), these are provided via the `sol` argument.

# Arguments
- `A::FEMatrix`: The finite element matrix to assemble into.
- `O::BilinearOperator`: The bilinear operator to assemble.
- `sol`: (Optional) Array of solution vectors or argument blocks required for nonlinear operators. Defaults to `nothing`.

# Keyword Arguments
- `assemble_matrix::Bool`: Whether to assemble the matrix (default: `true`).
- `kwargs...`: Additional keyword arguments passed to the assembly routines (e.g., for setting the time, verbosity, etc.).

"""
function assemble!(A::FEMatrix, O::BilinearOperator{Tv, UT}, sol = nothing; assemble_matrix = true, kwargs...) where {Tv, UT}
    if !assemble_matrix
        return nothing
    end
    ind_test = O.u_test
    ind_ansatz = O.u_ansatz
    ind_args = O.u_args
    return if length(O.u_args) > 0
        build_assembler!(A.entries, O, [A[j, j] for j in ind_test], [A[j, j] for j in ind_ansatz], [sol[j] for j in ind_args]; kwargs...)
        O.assembler(A.entries, nothing, [sol[j] for j in ind_args])
    else
        build_assembler!(A.entries, O, [A[j, j] for j in ind_test], [A[j, j] for j in ind_ansatz]; kwargs...)
        O.assembler(A.entries, nothing)
    end
end


function assemble!(A, b, sol, O::BilinearOperatorFromMatrix{UT, MT}, SC::SolverConfiguration; assemble_matrix = true, kwargs...) where {UT, MT}
    if !assemble_matrix
        return nothing
    end
    if UT <: Integer
        ind_test = O.u_test
        ind_ansatz = O.u_ansatz
    elseif UT <: Unknown
        ind_test = [get_unknown_id(SC, u) for u in O.u_test]
        ind_ansatz = [get_unknown_id(SC, u) for u in O.u_ansatz]
    end
    return if MT <: FEMatrix
        for (j, ij) in enumerate(ind_test), (k, ik) in enumerate(ind_ansatz)
            addblock!(A[j, k], O.A[ij, ik]; factor = O.parameters[:factor])
            if O.parameters[:transposed_copy] != 0
                addblock!(A[k, j], O.A[ij, ik]; transpose = true, factor = O.parameters[:factor] * O.parameters[:transposed_copy])
            end
        end
    else
        @assert length(ind_test) == 1
        @assert length(ind_ansatz) == 1
        addblock!(A[ind_test[1], ind_ansatz[1]], O.A; factor = O.parameters[:factor])
        if O.parameters[:transposed_copy] != 0
            addblock!(A[ind_ansatz[1], ind_test[1]], O.A; transpose = true, factor = O.parameters[:factor] * O.parameters[:transposed_copy])
        end
    end
end
