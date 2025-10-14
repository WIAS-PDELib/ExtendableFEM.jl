struct CoupledDofsRestriction{TM} <: AbstractRestriction
    coupling_matrices::Vector{TM}
    parameters::Dict{Symbol, Any}
end


"""
    CoupledDofsRestriction(matrix::AbstractMatrix)

    Creates a restriction that couples dofs together.

    The coupling is stored in a CSC Matrix `matrix`, s.t.,

    dofᵢ = Σⱼ Aⱼᵢ dofⱼ (col-wise)

    The matrix can be obtained from, e.g., `get_periodic_coupling_matrix`.
"""
function CoupledDofsRestriction(matrix::AbstractMatrix)
    return CoupledDofsRestriction([matrix], Dict{Symbol, Any}(:name => "CoupledDofsRestriction"))
end


"""
    CoupledDofsRestriction(matrices::Vector{AM}) where {AM <: AbstractMatrix}

    Creates a `CoupledDofsRestriction` from multiple given coupling matrices.
"""
function CoupledDofsRestriction(matrices::Vector{AM}) where {AM <: AbstractMatrix}
    return CoupledDofsRestriction(matrices, Dict{Symbol, Any}(:name => "CoupledDofsRestriction"))
end


function assemble!(R::CoupledDofsRestriction, sol, SC; kwargs...)

    # extract all col indices and remove duplicates
    # subtract diagonal and shrink matrix to non-empty cols
    Bs = [ (matrix - LinearAlgebra.I)[:, unique(findnz(matrix)[2])] for matrix in R.coupling_matrices ]

    # combine all into one matrix
    B = hcat(Bs...)

    # eliminate redundant cols by QR:
    qr_result = qr(B)

    # pick minimal number of cols that are rank preserving
    cols_of_interest = qr_result.pcol[1:rank(qr_result)]
    B = B[:, cols_of_interest]

    R.parameters[:matrix] = B
    R.parameters[:rhs] = Zeros(size(B, 2))

    # fixed dofs are all active rows of B
    R.parameters[:fixed_dofs] = unique(findnz(B)[1])

    return nothing
end
