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
    return CoupledDofsRestriction(
        [matrix],
        Dict{Symbol, Any}(
            :name => "CoupledDofsRestriction",
            :reduce_col_space => false
        )
    )
end


"""
    CoupledDofsRestriction(matrices::Vector{AM}) where {AM <: AbstractMatrix}

    Creates a `CoupledDofsRestriction` from multiple given coupling matrices.

    By default, the column space of the matrices is reduced to be of full rank.
    This can toggled by the `:reduce_col_space` parameter.
"""
function CoupledDofsRestriction(matrices::Vector{AM}) where {AM <: AbstractMatrix}
    return CoupledDofsRestriction(
        matrices,
        Dict{Symbol, Any}(
            :name => "CoupledDofsRestriction",
            :reduce_col_space => true
        )
    )
end


"""
    CoupledDofsRestriction(source::Ti, target::Ti, axis::Ti) where Ti

    Create an incomplete CoupledDofsRestriction by only providing abstract information:
    - unknown: the unknown to be coupled
    - source_region: source boundary region
    - target_region: target boundary region

    This constructor assumes that
    - the source and target boundary regions form two parallel hyperplanes
    - the coupling is performed orthogonally between the given hyperplanes

    The concrete coupling matrix will be computed in the solver, when the grid geometry and the FES is known.
"""
function CoupledDofsRestriction(unknown::Unknown, source_region::Ti, target_region::Ti; kwargs...) where {Ti}
    return CoupledDofsRestriction(
        [nothing],
        Dict{Symbol, Any}(
            :name => "CoupledDofsRestriction",
            :reduce_col_space => false,
            :unknown => unknown,
            :source_region => source_region,
            :target_region => target_region,
            :kwargs => kwargs
        )
    )
end


function assemble!(R::CoupledDofsRestriction, sol, SC; kwargs...)

    # remember original R
    R_orig = R

    # compute the coupling matrix first, if necessary
    if R.coupling_matrices[begin] === nothing
        FES = sol[R.parameters[:unknown]].FES
        source_region = R.parameters[:source_region]
        target_region = R.parameters[:target_region]
        R_kwargs = R.parameters[:kwargs]

        grid = FES.dofgrid

        # at first, we select one face on source/target region
        source_bface = findfirst(==(source_region), grid[BFaceRegions])
        target_bface = findfirst(==(target_region), grid[BFaceRegions])

        # the normal (it should really be opposite to the normal on the other side)
        normal = grid[BFaceNormals][:, source_bface]
        @assert normal ≈ -grid[BFaceNormals][:, target_bface]

        # pick a coordinate on each boundary region
        source_coord = grid[Coordinates][:, grid[BFaceNodes][1, source_bface]]
        target_coord = grid[Coordinates][:, grid[BFaceNodes][1, target_bface]]

        # compute the sum of scalar product with the normals (reflection point)
        γ = source_coord'normal + target_coord'normal

        function give_opposite!(y, x)
            σ = 2.0 * normal'x
            @. y = x + (γ - σ) * normal # then x ⇔ y are opposite along the normal vector
            return nothing
        end

        coupling_matrix = get_periodic_coupling_matrix(FES, source_region, target_region, give_opposite!; R_kwargs...)

        # replace R (in this scope)
        R = CoupledDofsRestriction(coupling_matrix)
    end


    # extract all col indices and remove duplicates
    # subtract diagonal and shrink matrix to non-empty cols
    Bs = [ (matrix - LinearAlgebra.I)[:, unique(findnz(matrix)[2])] for matrix in R.coupling_matrices ]

    # combine all into one matrix
    B = hcat(Bs...)

    if R.parameters[:reduce_col_space]
        # eliminate redundant cols by QR:
        qr_result = qr(B)

        # pick minimal number of cols that are rank preserving
        cols_of_interest = qr_result.pcol[1:rank(qr_result)]
        B = B[:, cols_of_interest]
    end

    R_orig.parameters[:matrix] = B
    R_orig.parameters[:rhs] = zeros(size(B, 2))
    R_orig.parameters[:multiplier] = zeros(size(B, 2))

    # fixed dofs are all active rows of B
    R_orig.parameters[:fixed_dofs] = unique(findnz(B)[1])

    return nothing
end
