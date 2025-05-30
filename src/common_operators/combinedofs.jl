###########################################
### COMBINE DOFS (e.g. for periodicity) ###
###########################################

mutable struct CombineDofs{UT, CT} <: AbstractOperator
    uX::UT                  # component nr for dofsX
    uY::UT                  # component nr for dofsY
    coupling_info::CT
    FESX::Any
    FESY::Any
    assembler::Any
    parameters::Dict{Symbol, Any}
end

default_combop_kwargs() = Dict{Symbol, Tuple{Any, String}}(
    :name => ("CombineDofs", "name for operator used in printouts"),
    :penalty => (1.0e30, "penalty for fixed degrees of freedom"),
    :verbosity => (0, "verbosity level"),
)

# informs solver in which blocks the operator assembles to
function dependencies_when_linearized(O::CombineDofs)
    return [O.uX, O.uY]
end


"""
````
function CombineDofs(uX, uY, dofsX, dofsY, factors; kwargs...)
````

When assembled, the dofsX of the unknown uX will be coupled
with the dofsY of uY, e.g., for periodic boundary conditions.
 

Keyword arguments:
$(_myprint(default_combop_kwargs()))

"""
function CombineDofs(uX, uY, dofsX, dofsY, factors = ones(Int, length(dofsX)); kwargs...)
    # build a sparse matrix from dofsX, dofsY, factors
    # we have to set the size (maximum of indices: matrix may be smaller than the system matrix, this is ok)
    # for sparse matrices, define a `combine` function for duplicates: pick the latest entry
    # convention: dofsX acts as column indices!
    coupling_matrix = sparse(dofsY, dofsX, factors, maximum(dofsY), maximum(dofsX), (a, b) -> b)
    return CombineDofs(uX, uY, coupling_matrix; kwargs...)
end

"""
````
function CombineDofs(uX, uY, coupling_matrix::AbstractMatrix; kwargs...)
````
Input:
    uX: FEVectorBlock of the "from" boundary region
    uY: FEVectorBlock of the "to" boundary region (usually the same as uX)
    coupling_matrix: coupling matrix computed from `get_periodic_coupling_matrix`

Keyword arguments:
$(_myprint(default_combop_kwargs()))
"""
function CombineDofs(uX, uY, coupling_matrix::AbstractMatrix; kwargs...)
    parameters = Dict{Symbol, Any}(k => v[1] for (k, v) in default_combop_kwargs())
    _update_params!(parameters, kwargs)
    return CombineDofs(uX, uY, coupling_matrix, nothing, nothing, nothing, parameters)
end

function apply_penalties!(A, b, sol, CD::CombineDofs{UT, CT}, SC::SolverConfiguration; assemble_matrix = true, assemble_rhs = true, kwargs...) where {UT, CT}
    if UT <: Integer
        ind = [CD.ux, CD.uY]
    elseif UT <: Unknown
        ind = [get_unknown_id(SC, CD.uX), get_unknown_id(SC, CD.uY)]
    end
    build_assembler!(CD, [sol[j] for j in ind])
    return CD.assembler(A.entries, b.entries, assemble_matrix, assemble_rhs)
end


function build_assembler!(CD::CombineDofs{UT, CT}, FE::Array{<:FEVectorBlock, 1}; time = 0.0) where {UT, CT <: AbstractMatrix}
    ## check if FES is the same as last time
    FESX, FESY = FE[1].FES, FE[2].FES
    if (CD.FESX != FESX) || (CD.FESY != FESY)
        coupling_matrix = CD.coupling_info
        offsetX = FE[1].offset
        offsetY = FE[2].offset
        if CD.parameters[:verbosity] > 0
            @info ".... coupling $(length(coupling_matrix.nzval)) dofs"
        end
        penalty = CD.parameters[:penalty]

        function assemble!(A::AbstractSparseArray{T}, b::AbstractVector{T}, assemble_matrix::Bool, assemble_rhs::Bool, kwargs...) where {T}

            timerOutput = TimerOutput()

            if assemble_matrix
                # go through each constrained dof and update the FE adjacency info
                # of the coupled dofs
                @timeit timerOutput "for loop 1" for dof_i in 1:size(coupling_matrix, 2)

                    # this col-view is efficient
                    @timeit timerOutput "coupling_i view" coupling_i = @views coupling_matrix[:, dof_i]

                    # do nothing if dof_k is not coupled to any constrained dof
                    @timeit timerOutput "nnz == 0" if nnz(coupling_i) == 0
                        continue
                    end

                    # write the FE adjacency of the constrained dofs into this row
                    sourcerow = dof_i + offsetX

                    # extract the constrained dofs and the weights
                    @timeit timerOutput "findnz" coupled_dofs_i, weights_i = findnz(coupling_i)

                    # parse through sourcerow and add the contents to the coupled dofs
                    @timeit timerOutput  "A[sourcerow, col]"  for col in 1:size(A, 2)
                        val = A[sourcerow, col]
                        if abs(val) > 1.0e-15
                            for (dof_k, weight_ik) in zip(coupled_dofs_i, weights_i)
                                targetrow = dof_k + offsetX
                                _addnz(A, targetrow, col, val, weight_ik)
                            end
                        end
                    end
                end

                # replace the geometric coupling rows based
                # on the original coupling matrix
                @timeit timerOutput "for loop 2" for dof_i in 1:size(coupling_matrix, 2)

                    coupling_i = coupling_matrix[:, dof_i]
                    # do nothing if no coupling for dof_i
                    if nnz(coupling_i) == 0
                        continue
                    end

                    # get the coupled dofs of dof_i and the corresponding weights
                    @timeit timerOutput "findnz"  coupled_dofs_i, weights_i = findnz(coupling_i)

                    sourcerow = dof_i + offsetX

                    # replace sourcerow with coupling linear combination
                    @timeit timerOutput "_addnz"  begin
                        _addnz(A, sourcerow, sourcerow, -1.0, penalty)
                        for (dof_j, weight_ij) in zip(coupled_dofs_i, weights_i)
                            # weights for ∑ⱼ wⱼdofⱼ - dofᵢ = 0
                            _addnz(A, sourcerow, dof_j + offsetY, weight_ij, penalty)
                        end
                    end
                end
                @timeit timerOutput "flush"  flush!(A)
            end

            if assemble_rhs

                for dof_i in 1:size(coupling_matrix, 2)
                    # this col-view is efficient
                    coupling_i = @views coupling_matrix[:, dof_i]
                    # do nothing if no coupling for dof_i
                    if nnz(coupling_i) == 0
                        continue
                    end

                    # get the coupled dofs of dof_i and the corresponding weights
                    coupled_dofs, weights = findnz(coupling_i)

                    # transfer all assembly information to dof_i
                    sourcerow = dof_i + offsetY
                    for (dof_j, weight) in zip(coupled_dofs, weights)
                        targetrow = dof_j + offsetY
                        b[targetrow] += weight * b[sourcerow]
                    end
                end


                # now set the rows of the constrained dofs to zero to enforce the linear combination
                for dof_i in 1:size(coupling_matrix, 2)
                    coupling_i = coupling_matrix[:, dof_i]
                    # do nothing if no coupling for dof_i
                    if nnz(coupling_i) == 0
                        continue
                    end

                    b[dof_i + offsetX] = 0.0

                end
            end

            TimerOutputs.complement!(timerOutput)
            show(timerOutput)

            return nothing
        end

        CD.assembler = assemble!
        CD.FESX = FESX
        CD.FESY = FESY
    end

    return nothing
end
