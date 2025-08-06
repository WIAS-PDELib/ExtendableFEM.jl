struct CoupledDofsRestriction{TM} <: AbstractRestriction
    coupling_matrix::TM
    parameters::Dict{Symbol, Any}
end


"""
    CoupledDofsRestriction(matrix::AbstractMatrix)

    Creates an restriction that couples dofs together.

    The coupling is stored in a CSC Matrix `matrix`, s.t.,

    dofᵢ = Σⱼ Aⱼᵢ dofⱼ (col-wise)

    The matrix can be obtained from, e.g., `get_periodic_coupling_matrix`.
"""
function CoupledDofsRestriction(matrix::AbstractMatrix)
    return CoupledDofsRestriction(matrix, Dict{Symbol, Any}(:name => "CoupledDofsRestriction"))
end


function assemble!(R::CoupledDofsRestriction, SC; kwargs...)

    # extract all col indices
    _, J, _ = findnz(R.coupling_matrix)

    # remove duplicates
    unique_cols = unique(J)

    # subtract diagonal and shrink matrix to non-empty cols
    B = (R.coupling_matrix - LinearAlgebra.I)[:, unique_cols]

    R.parameters[:matrix] = B
    R.parameters[:rhs] = Zeros(length(unique_cols))

    return nothing
end
