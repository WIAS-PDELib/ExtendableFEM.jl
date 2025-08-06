struct CoupledDofsRestriction{TM} <: AbstractRestriction
    coupling_matrix::TM
    parameters::Dict{Symbol, Any}
end

function CoupledDofsRestriction(matrix::AbstractMatrix)
    return CoupledDofsRestriction(matrix, Dict{Symbol, Any}(:name => "CoupledDofsRestriction"))
end


function assemble!(R::CoupledDofsRestriction, A, b, sol, SC; kwargs...)

    # extract all col indices
    _, J, _ = findnz(R.coupling_matrix)

    # remove duplicates
    unique_cols = unique(J)

    # subtract diagonal and shrink matrix to non-empty cols
    B = (R.coupling_matrix - LinearAlgebra.I)[:, unique_cols]

    R.parameters[:matrix] = B
    return R.parameters[:rhs] = Zeros(length(unique_cols))

    return nothing
end
