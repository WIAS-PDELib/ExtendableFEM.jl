###########################################
### COMBINE DOFS (e.g. for periodicity) ###
###########################################

mutable struct CombineDofs{UT, CT} <: AbstractOperator
    uX::UT                  # component nr for dofsX
    uY::UT                  # component nr for dofsY
    coupling_info::CT
    #dofsX::XT
    #dofsY::YT
    #factors::FT
    FESX::Any
    FESY::Any
    assembler::Any
    parameters::Dict{Symbol, Any}
end

default_combop_kwargs() = Dict{Symbol, Tuple{Any, String}}(
    :penalty => (1.0e30, "penalty for fixed degrees of freedom"),
    :verbosity => (0, "verbosity level"),
)

# informs solver in which blocks the operator assembles to
function dependencies_when_linearized(O::CombineDofs)
    return [O.uX, O.uY]
end

function fixed_dofs(O::CombineDofs{UT, CT}) where {UT, CT <: Tuple}
    ## assembles operator to full matrix A and b
    return O.coupling_info[2]
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
function CombineDofs(uX, uY, dofsX, dofsY, factors = ones(Int, length(X)); kwargs...)
    parameters = Dict{Symbol, Any}(k => v[1] for (k, v) in default_combop_kwargs())
    _update_params!(parameters, kwargs)
    @assert length(dofsX) == length(dofsY)
    coupling_info = (dofsX, dofsY, factors)
    return CombineDofs(uX, uY, coupling_info, nothing, nothing, nothing, parameters)
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
    @assert size(coupling_matrix, 1) == size(coupling_matrix, 2)
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

function build_assembler!(CD::CombineDofs{UT, TupleType}, FE::Array{<:FEVectorBlock, 1}; time = 0.0) where {UT, TupleType <: Tuple}
    ## check if FES is the same as last time
    FESX, FESY = FE[1].FES, FE[2].FES
    if (CD.FESX != FESX) || (CD.FESY != FESY)
        dofsX = CD.coupling_info[1]
        dofsY = CD.coupling_info[2]
        factors = CD.coupling_info[3]
        offsetX = FE[1].offset
        offsetY = FE[2].offset
        if CD.parameters[:verbosity] > 0
            @info ".... combining $(length(dofsX)) dofs"
        end
        function assemble(A::AbstractSparseArray{T}, b::AbstractVector{T}, assemble_matrix::Bool, assemble_rhs::Bool, kwargs...) where {T}
            if assemble_matrix
                targetrow::Int = 0
                sourcerow::Int = 0
                targetcol::Int = 0
                sourcecol::Int = 0
                val::Float64 = 0
                ncols::Int = size(A, 2)
                for gdof in eachindex(dofsX)
                    # copy source row (for dofY) to target row (for dofX)
                    targetrow = dofsX[gdof] + offsetX
                    sourcerow = offsetY + dofsY[gdof]
                    for sourcecol in 1:ncols
                        targetcol = sourcecol - offsetY + offsetX
                        val = A[sourcerow, sourcecol]
                        _addnz(A, targetrow, targetcol, factors[gdof] * val, 1)
                        A[sourcerow, sourcecol] = 0
                    end

                    # replace source row (of dofY) with equation for coupling the two dofs
                    sourcecol = dofsY[gdof] + offsetY
                    targetcol = dofsX[gdof] + offsetX
                    sourcerow = offsetY + dofsY[gdof]
                    _addnz(A, sourcerow, targetcol, 1, 1)
                    _addnz(A, sourcerow, sourcecol, -factors[gdof], 1)
                end
                flush!(A)
            end
            return if assemble_rhs
                for gdof in 1:length(dofsX)
                    sourcerow = offsetY + dofsY[gdof]
                    targetrow = offsetX + dofsX[gdof]
                    b[targetrow] += b[sourcerow]
                    b[sourcerow] = 0
                end
            end
        end
        CD.assembler = assemble
        CD.FESX = FESX
        CD.FESY = FESY
    end

    return nothing
end


function build_assembler!(CD::CombineDofs{UT, MatrixType}, FE::Array{<:FEVectorBlock, 1}; time = 0.0) where {UT, MatrixType <: AbstractMatrix}
    ## check if FES is the same as last time
    FESX, FESY = FE[1].FES, FE[2].FES
    if (CD.FESX != FESX) || (CD.FESY != FESY)
        coupling_matrix = CD.coupling_info
        offsetX = FE[1].offset
        offsetY = FE[2].offset
        if CD.parameters[:verbosity] > 0
            @info ".... coupling $(length(coupling_matrix.nzval)) dofs"
        end
        function assemble(A::AbstractSparseArray{T}, b::AbstractVector{T}, assemble_matrix::Bool, assemble_rhs::Bool, kwargs...) where {T}
            if assemble_matrix
                for dof_i in 1:size(coupling_matrix, 2)
                    # this col-view is efficient
                    coupling_i = @views coupling_matrix[:, dof_i]
                    # do nothing if no coupling for dof_i
                    if nnz(coupling_i) == 0
                        continue
                    end

                    # get the coupled dofs of dof_i and the corresponding weights
                    coupled_dofs, weights = findnz(coupling_i)

                    # transfer all assembly information of dof_i to the
                    # coupled dofs: add the corresponding matrix row
                    sourcerow = dof_i + offsetY
                    for sourcecol in 1:size(A, 2)
                        val = A[sourcerow, sourcecol]
                        for dof_j in coupled_dofs
                            targetrow = dof_j + offsetX
                            targetcol = sourcecol
                            _addnz(A, targetrow, targetcol, val, 1)
                        end
                        # eliminate the sourcerow
                        A[sourcerow, sourcecol] = 0
                    end


                    # replace sourcerow with coupling linear combination
                    _addnz(A, sourcerow, sourcerow, -1.0, 1)
                    for (dof_j, weight) in zip(coupled_dofs, weights)
                        # set negative weight for dofᵢ - ∑ⱼ wⱼdofⱼ = 0
                        _addnz(A, sourcerow, dof_j + offsetY, weight, 1)
                    end

                    flush!(A)
                end
            end

            if assemble_rhs
                for dof_i in 1:length(b)
                    # this col-view is efficient
                    coupling_i = @views coupling_matrix[:, dof_i]
                    # do nothing if no coupling for dof_i
                    if nnz(coupling_i) == 0
                        continue
                    end

                    coupled_dofs, weights = findnz(coupling_i)
                    sourcerow = dof_i + offsetY

                    for dof_j in coupled_dofs
                        targetrow = dof_j + offsetX
                        b[targetrow] += b[sourcerow]
                    end

                    b[sourcerow] = 0.0
                end
            end

            return nothing
        end

        CD.assembler = assemble
        CD.FESX = FESX
        CD.FESY = FESY
    end

    return nothing
end
