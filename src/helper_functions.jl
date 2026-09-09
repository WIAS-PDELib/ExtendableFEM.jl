"""
````
function get_periodic_coupling_info(FES, xgrid, b1, b2, is_opposite::Function; factor_vectordofs = "auto")
````

computes the dofs that have to be coupled for periodic boundary conditions on the given xgrid for boundary regions b1, b2.
The is_opposite function evaluates if two provided face midpoints are on opposite sides to each other (the mesh xgrid should be appropriate).
For vector-valued FETypes the user can provide factor_vectordofs to incorporate a sign change if needed.
This is automatically done for all Hdiv-conforming elements and (for the normal-weighted face bubbles of) the Bernardi-Raugel element H1BR. 

"""
function get_periodic_coupling_info(
        FES::FESpace,
        xgrid::ExtendableGrid,
        b1,
        b2,
        is_opposite::Function;
        factor_vectordofs = "auto",
        factor_components = "auto"
    )

    FEType = eltype(FES)
    ncomponents = get_ncomponents(FEType)
    if factor_vectordofs == "auto"
        if FEType <: AbstractHdivFiniteElement || FEType <: H1BR
            factor_vectordofs = -1
        else
            factor_vectordofs = 1
        end
    end
    if factor_components == "auto"
        factor_components = ones(Int, ncomponents)
    end


    @assert FEType <: AbstractH1FiniteElement "not yet working for non H1-conforming elements"
    xBFaceRegions = xgrid[BFaceRegions]
    xBFaceNodes = xgrid[BFaceNodes]
    xBFaceFaces = xgrid[BFaceFaces]
    xCoordinates = xgrid[Coordinates]
    nbfaces = size(xBFaceNodes, 2)
    nnodes = num_nodes(xgrid)
    nnodes4bface = size(xBFaceNodes, 1)
    EG = xgrid[UniqueBFaceGeometries][1]
    xdim = size(xCoordinates, 1)
    nedges4bface = xdim == 3 ? num_faces(EG) : 0
    xBFaceMidPoints = zeros(Float64, xdim, nbfaces)
    for bface in 1:nbfaces, j in 1:xdim, bn in 1:nnodes4bface
        xBFaceMidPoints[j, bface] += xCoordinates[j, xBFaceNodes[bn, bface]] / nnodes4bface
    end
    if xdim == 3
        xEdgeMidPoint = zeros(Float64, xdim)
        xEdgeMidPoint2 = zeros(Float64, xdim)
        xEdgeNodes = xgrid[EdgeNodes]
        xFaceEdges = xgrid[FaceEdges]
        if FEType <: H1P1
            nedgedofs = 0
        elseif FEType <: H1P2
            nedgedofs = 1
        elseif FEType <: H1P3
            nedgedofs = 2
        else
            @warn "get_periodic_coupling_info not yet working for non H1-conforming elements"
        end
    end
    @assert FEType <: AbstractH1FiniteElement "get_periodic_coupling_info not yet working for non H1-conforming elements"
    xBFaceDofs = FES[BFaceDofs]
    dofsX, dofsY, factors = Int[], Int[], Int[]
    counterface = 0
    nfb = 0
    partners = zeros(Int, xdim)
    coffsets = ExtendableFEMBase.get_local_coffsets(FEType, ON_BFACES, EG)
    nedgedofs = 0

    for bface in 1:nbfaces
        counterface = 0
        nfb = num_targets(xBFaceDofs, bface)
        if xBFaceRegions[bface] == b1
            for bface2 in 1:nbfaces
                if xBFaceRegions[bface2] == b2
                    if is_opposite(view(xBFaceMidPoints, :, bface), view(xBFaceMidPoints, :, bface2))
                        counterface = bface2
                        break
                    end
                end
            end
        end
        if counterface > 0

            # couple first two node dofs in opposite order due to orientation
            for c in 1:ncomponents
                if factor_components[c] == 0
                    continue
                end
                nfbc = coffsets[c + 1] - coffsets[c] # total dof count for this component

                # couple nodes
                for nb in 1:nnodes4bface
                    ## find node partner on other side that evaluates true in is_ooposite function
                    for nc in 1:nnodes4bface
                        if is_opposite(view(xCoordinates, :, xBFaceDofs[nb, bface]), view(xCoordinates, :, xBFaceDofs[nc, counterface]))
                            partners[nb] = nc
                            break
                        end
                    end
                    ## couple node dofs (to be skipped for e.g. Hdiv, Hcurl elements)
                    push!(dofsX, xBFaceDofs[coffsets[c] + nb, bface])
                    push!(dofsY, xBFaceDofs[coffsets[c] + partners[nb], counterface])
                end
                # @info "matching face $bface (nodes = $(xBFaceNodes[:,bface]), dofs = $(xBFaceDofs[:,bface])) with face $counterface (nodes = $(xBFaceNodes[:,counterface]), dofs = $(xBFaceDofs[:,counterface])) with partner node order $partners"

                ## couple edges
                if nedges4bface > 0 && FEType <: H1P2 || FEType <: H1P3
                    # todo: for H1P3 edge orientation place a role !!!
                    for nb in 1:nedges4bface
                        fill!(xEdgeMidPoint, 0)
                        for j in 1:xdim, k in 1:2
                            xEdgeMidPoint[j] += xCoordinates[j, xEdgeNodes[k, xFaceEdges[nb, xBFaceFaces[bface]]]] / 2
                        end
                        ## find edge partner on other side that evaluates true at edge midpoint in is_opposite function
                        for nc in 1:nnodes4bface
                            fill!(xEdgeMidPoint2, 0)
                            for j in 1:xdim, k in 1:2
                                xEdgeMidPoint2[j] += xCoordinates[j, xEdgeNodes[k, xFaceEdges[nc, xBFaceFaces[counterface]]]] / 2
                            end
                            if is_opposite(xEdgeMidPoint, xEdgeMidPoint2)
                                partners[nb] = nc
                                break
                            end
                        end

                        ## couple edge dofs (local orientation information is needed for more than one dof on each edge !!! )
                        for k in 1:nedgedofs
                            push!(dofsX, xBFaceDofs[coffsets[c] + nnodes4bface + nb + (k - 1) * nedgedofs, bface])
                            push!(dofsY, xBFaceDofs[coffsets[c] + nnodes4bface + partners[nb] + (k - 1) * nedgedofs, counterface])
                        end
                    end
                end

                ## couple face dofs (interior dofs of bface)
                for nb in 1:(nfbc - nnodes4bface - nedges4bface * nedgedofs)
                    push!(dofsX, xBFaceDofs[coffsets[c] + nnodes4bface + nedges4bface * nedgedofs + nb, bface])
                    push!(dofsY, xBFaceDofs[coffsets[c] + nnodes4bface + nfbc - nnodes4bface - nedges4bface * nedgedofs + 1 - nb, counterface]) # couple face dofs in opposite order due to orientation (works in 2D at least)
                end
                append!(factors, ones(nfbc) * factor_components[c])
            end

            ## couple remaining dofs (should be vector dofs)
            for dof in (coffsets[end] + 1):nfb
                push!(dofsX, xBFaceDofs[dof, bface])
                push!(dofsY, xBFaceDofs[nfb - coffsets[end] + dof - 1, counterface]) # couple face dofs in opposite order due to orientation (works in 2D at least, e.g. for Bernardi--Raugel)
                push!(factors, factor_vectordofs)
            end
        end
    end

    return dofsX, dofsY, factors
end


# compact variant of lazy_interpolate! specialized on ON_FACES interpolations
function interpolate_on_boundaryfaces(
        source::FEVector{Tv, TvG, TiG},
        xgrid::ExtendableGrid{TvG, TiG},
        source_target_transform!,
        post_mutation!,
        start_cell::Int = 1, # TODO we interpolate on the "b_from" side: a proper start cell should be given
        eps = 1.0e-13,
        kwargs...,
    ) where {Tv, TvG, TiG}

    # wrap point evaluation into function that is put into normal interpolate!
    xdim::Int = size(xgrid[Coordinates], 1)
    PE = PointEvaluator(post_mutation!, [(1, Identity)], source)
    xref = zeros(TvG, xdim)
    x_source = zeros(TvG, xdim)
    CF::ExtendableGrids.CellFinder{TvG, TiG} = ExtendableGrids.CellFinder(xgrid)
    last_cell = zeros(Int, 1)
    last_cell[1] = start_cell

    function __setstartcell(new)
        return last_cell[1] = Int(new)
    end

    function __eval_point(result, qpinfo)
        source_target_transform!(x_source, qpinfo.x)

        cell = ExtendableGrids.gFindLocal!(xref, CF, x_source; icellstart = last_cell[1], eps)
        if cell == 0
            @error "boundary coordinate $(qpinfo.x) opposite to $x_source could not be found in the grid"
        else
            evaluate_bary!(result, PE, xref, cell)
            last_cell[1] = cell
        end
        return nothing
    end

    return __eval_point, __setstartcell
end

function get_periodic_coupling_matrix(
        FES::FESpace{Tv},
        xgrid::ExtendableGrid{TvG, TiG},
        b_from,
        b_to,
        source_target_transform!::Function;
        post_mutation! = ExtendableFEMBase.standard_kernel,
        kwargs...
    ) where {Tv, TvG, TiG}
    @warn "get_periodic_coupling_matrix with grid argument is deprecated"
    return _get_periodic_coupling_matrix(FES, xgrid, b_from, b_to, source_target_transform!, post_mutation!; kwargs...)
end

# merge matrix B into A, overriding the entries of A if an entry is both present in A and B
function merge!(A::ExtendableSparseMatrix, B::ExtendableSparseMatrix)
    rows, cols, values = findnz(B)
    for (row, col, value) in zip(rows, cols, values)
        A[row, col] = value
    end
    return nothing
end

function _get_periodic_coupling_matrix(
        FES::FESpace{Tv},
        xgrid::ExtendableGrid{TvG, TiG},
        b_from,
        b_to,
        source_target_transform!::Function,
        post_mutation!::Function;
        mask = :auto,
        sparsity_tol = 1.0e-12,
        parallel = false,
        threads = Threads.nthreads()
    ) where {Tv, TvG, TiG}

    nthr = parallel ? threads : 1
    @info "Computing periodic coupling matrix with $nthr thread(s). This may take a while."

    if typeof(b_from) <: Int
        b_from = [b_from]
    end

    if typeof(b_to) <: Int
        b_to = [b_to]
    end

    # total number of grid boundary faces
    boundary_nodes = xgrid[BFaceNodes]
    n_boundary_faces = size(boundary_nodes, 2)

    # corresponding boundary regions
    boundary_regions = xgrid[BFaceRegions]

    # FE basis dofs on each boundary face
    dofs_on_boundary = FES[BFaceDofs]

    # find all faces in b_to
    bfaces_in_b_to = TiG[]
    bfaces_in_b_from = TiG[]
    for (i, region) in enumerate(boundary_regions)
        if region in b_to
            push!(bfaces_in_b_to, i)
        elseif region in b_from
            push!(bfaces_in_b_from, i)
        end
    end

    # offset of the individual components of the FES
    ncomponents = get_ncomponents(get_FEType(FES))
    coffset = FES.coffset

    # fill component mask if not done before
    if mask == :auto
        mask = ones(ncomponents)
    else
        @assert length(mask) == ncomponents "component mask has to match number of components"
    end

    # do the intervals a=[a1,a2] and b=[b1,b2] overlap?
    safety = 1.0e-12 # we would be sad if we miss an overlap due to rounding errors
    do_intervals_overlap(a, b) = a[1] ≤ b[2] + safety && b[1] ≤ a[2] + safety

    # do the boxes 𝑓 and 𝑔 overlap?
    # we provide the Vector of the coordinate intervals
    function do_boxes_overlap(box_f::AbstractVector, box_g::AbstractVector)
        for i in eachindex(box_f)
            if !do_intervals_overlap(box_f[i], box_g[i])
                return false
            end
        end

        # all coordinates overlap
        return true
    end

    dummy = zeros(TvG, size(xgrid[Coordinates], 1))

    # transform a face to the target side using the source_target_transform! function
    # Warning: this overwrites the face
    function transfer_face!(face::AbstractMatrix)
        for i in axes(face, 2)
            @views coord = face[:, i]
            source_target_transform!(dummy, coord)
            @views face[:, i] .= dummy
        end
        return
    end

    # precompute approximate search region for each boundary face in b_from
    searchareas = ExtendableGrids.VariableTargetAdjacency(TiG)
    coords = xgrid[Coordinates]
    bfacenodes = xgrid[BFaceNodes]
    coords_to = coords[:, bfacenodes[:, 1]]
    nodes_per_faces = size(coords_to, 2)
    dim = size(coords_to, 1)
    box_to = @MArray [Float64[0, 0], Float64[0, 0], Float64[0, 0]]
    for bface_to in bfaces_in_b_to
        for j in 1:nodes_per_faces, k in 1:dim
            coords_to[k, j] = coords[k, bfacenodes[j, bface_to]]
        end

        # transfer the coords_to to the other side
        transfer_face!(coords_to)

        # get the extrama in each component ( = bounding box of the face)
        for k in 1:dim
            box_to[k][1] = minimum(view(coords_to, k, :))
            box_to[k][2] = maximum(view(coords_to, k, :))
        end

        function inner_loop(bfaces_chunk)
            # some data
            local coords_from = coords[:, bfacenodes[:, 1]]
            local box_from = @MArray [Float64[0, 0], Float64[0, 0], Float64[0, 0]]
            local bfaces_from = Int[]

            for bface_from in bfaces_chunk
                for j in 1:nodes_per_faces, k in 1:dim
                    coords_from[k, j] = coords[k, bfacenodes[j, bface_from]]
                end
                for k in 1:dim
                    box_from[k][1] = minimum(view(coords_from, k, :))
                    box_from[k][2] = maximum(view(coords_from, k, :))
                end

                if do_boxes_overlap(box_to, box_from)
                    push!(bfaces_from, bface_from)
                end
            end

            return bfaces_from
        end

        if parallel && nthr > 1
            # create chunks to split this range for the threads
            bfaces_chunks = chunks(bfaces_in_b_from, n = nthr)

            tasks = map(bfaces_chunks) do bfaces_chunk
                Threads.@spawn inner_loop(bfaces_chunk)
            end

            # put all results together
            bfaces_from = vcat(fetch.(tasks)...)
        else
            bfaces_from = inner_loop(bfaces_in_b_from)
        end

        append!(searchareas, bfaces_from)
    end

    # flip the adjacency: in the following we need search areas for each "from" face
    searchareas = ExtendableGrids.atranspose(searchareas)

    # we are only interest in global bface numbers on the "from" boundary
    bfaces_of_interest = filter(bface -> boundary_regions[bface] in b_from, 1:n_boundary_faces)

    # throw error if no search area had been found for a bface
    for source in bfaces_of_interest
        if num_targets(searchareas, source) == 0
            throw("bface $source has no valid search area on the target side of the grid. Double check the provided from/to regions and your source_target_transform! function")
        end
    end

    # loop over boundary face indices in a chunk: we need this index for dofs_on_boundary
    function compute_chunk_result(chunk)

        # prepare data for this chunk
        # we need our own copy of the FE Space to avoid data race in the pre-computed interpolators
        our_FES = deepcopy(FES)
        local fe_vector = FEVector(our_FES)
        # to be sure
        fill!(fe_vector.entries, 0.0)

        local entries = SparseVector{Float64, Int64}(our_FES.ndofs, Int64[], Float64[])
        local fe_vector_target = FEVector(our_FES; entries)

        local n = length(fe_vector.entries)
        local result = ExtendableSparseMatrix(n, n)

        local eval_point, _ = interpolate_on_boundaryfaces(fe_vector, xgrid, source_target_transform!, post_mutation!)

        for boundary_face in chunk

            local local_dofs = @views dofs_on_boundary[:, boundary_face]
            for local_dof in local_dofs
                # compute number of component
                if mask[1 + ((local_dof - 1) ÷ coffset)] == 0.0
                    continue
                end

                # reset
                empty!(fe_vector_target.entries.nzind)
                empty!(fe_vector_target.entries.nzval)

                # activate one entry
                fe_vector.entries[local_dof] = 1.0

                interpolate!(
                    fe_vector_target[1],
                    ON_BFACES, eval_point,
                    items = @views bfaces_in_b_to[searchareas[:, boundary_face]]
                )

                # deactivate entry
                fe_vector.entries[local_dof] = 0.0

                # set entries
                for (i, target_entry) in zip(findnz(fe_vector_target.entries)...)
                    if abs(target_entry) > sparsity_tol
                        result[local_dof, i] = target_entry
                    end
                end
            end
        end

        return result
    end

    if parallel && nthr > 1
        # we loop ober the n_boundary_faces in parallel:
        # create chunks to split this range for the threads
        bface_chunks = chunks(bfaces_of_interest, n = nthr)

        # show start all chunks in parallel
        tasks = map(bface_chunks) do chunk
            Threads.@spawn compute_chunk_result(chunk)
        end

        # @info "done..."
        # wait for all chunks to finish and get results
        results = fetch.(tasks)

        # merge all matrices
        result = results[begin]
        for res_i in results[(begin + 1):end]
            merge!(result, res_i)
        end
    else
        result = compute_chunk_result(bfaces_of_interest)
    end

    sp_result = sparse(result)

    # strange if nothing is coupled
    if nnz(sp_result) == 0
        @warn "no coupling found. Are the grid boundary regions and the source_target_transform! method correct?"
    end

    return sp_result
end

"""
    get_periodic_coupling_matrix(
        FES::FESpace,
        b_from,
        b_to,
        source_target_transform!::Function;
        post_mutation! = ExtendableFEMBase.standard_kernel,
        mask = :auto,
        sparsity_tol = 1.0e-12
    )

Compute a coupling information for each dof on one boundary as a linear combination of dofs on another boundary

Input:
 - FES: FE space to be coupled (on its dofgrid)
 - b_from: boundary region(s) of the grid which dofs should be replaced in terms of dofs on b_to
 - b_to: boundary region(s) of the grid with dofs to replace the dofs in b_from
  - source_target_transform! Function `source_target_transform!(y, x)` that maps a point `x ∈ b_from` to the corresponding point `y` on the target boundary
 - post_mutation!: optional post-transformation applied after interpolation (e.g. to flip velocity components for Stokes problems)
 - mask: (optional) vector of masking components
 - sparsity_tol: threshold for treating an interpolated value as zero

source_target_transform!(x, y) has to be defined in a way that for each x ∈ b_from the resulting y is on the target boundary.
For each x in the grid, the resulting y has to be in the grid, too: incorporate some mirroring of the coordinates.
Example: If b_from is at x[1] = 0 and the target boundary is at y[1] = 1, then source_target_transform!(x, y) = y .= [ 1-x[1], x[2] ]

The return value is a (𝑛 × 𝑛) sparse matrix 𝐴 (𝑛 is the total number of dofs) containing the periodic coupling information.
The relation ship between the degrees of freedom is  dofᵢ = ∑ⱼ Aⱼᵢ ⋅ dofⱼ.
It is guaranteed that
    i)  Aᵢⱼ=0 if dofᵢ is 𝑛𝑜𝑡 on the boundary b_from.
    ii) Aᵢⱼ=0 if the opposite of dofᵢ is not in the same grid cell as dofⱼ.
"""
function get_periodic_coupling_matrix(
        FES,
        b_from,
        b_to,
        source_target_transform!;
        post_mutation! = ExtendableFEMBase.standard_kernel,
        kwargs...
    )
    return _get_periodic_coupling_matrix(FES, FES.dofgrid, b_from, b_to, source_target_transform!, post_mutation!; kwargs...)
end


## determines a common assembly grid for the given arrays of finite element spaces
function determine_assembly_grid(FES_test, FES_ansatz = [], FES_args = [])
    xgrid = FES_test[1].xgrid
    dofgrid = FES_test[1].dofgrid
    all_same_xgrid = true
    all_same_dofgrid = true
    for j in 2:length(FES_test)
        if xgrid !== FES_test[j].xgrid
            all_same_xgrid = false
        end
        if dofgrid !== FES_test[j].dofgrid
            all_same_dofgrid = false
        end
    end
    for j in 1:length(FES_ansatz)
        if xgrid !== FES_ansatz[j].xgrid
            all_same_xgrid = false
        end
        if dofgrid !== FES_ansatz[j].dofgrid
            all_same_dofgrid = false
        end
    end
    for j in 1:length(FES_args)
        if xgrid !== FES_args[j].xgrid
            all_same_xgrid = false
        end
        if dofgrid !== FES_args[j].dofgrid
            all_same_dofgrid = false
        end
    end
    if all_same_dofgrid
        return dofgrid
    elseif all_same_xgrid
        return xgrid
    else
        @warn "detected non-matching grids for involved finite element spaces, trying assembly on grid of first testfunction argument"
        return xgrid
    end
    return xgrid
end

## gets the dofmap for the FESpace FES fr the assemblygrid xgrid and the assembly type AT
function get_dofmap(FES, xgrid, AT)
    DM = Dofmap4AssemblyType(AT)
    if FES.dofgrid !== xgrid && FES.xgrid !== xgrid
        @warn "warning assembly grid does neither match FES dofgrid or parent grid!"
        return FES[DM]
    end
    FES[DM]
    return FES.dofgrid === xgrid ? FES[DM] : FES[ParentDofmap4Dofmap(DM)]
end


# """
# ````
# function tensor_view(input, i, rank, dim)
# ````

# Returns a view of input[i] and following entries
# reshaped as a tensor of rank `rank`.
# The parameter `dim` specifies the size of a tensor in each direction,
# e.g. a 1-Tensor (Vector) of length(dim) or a dim x dim 2-Tensor (matrix).
# As an example `tensor_view(v,5,2,5)` returns a view of `v(5:29)`
# as a 5x5 matrix.

# """
# function tensor_view(input, i::Int, rank::Int, dim::Int)
#     if rank == 0
#         return view(input, i:i)
#     elseif rank == 1
#         return view(input, i:i+dim-1)
#     elseif rank == 2
#         return reshape(view(input, i:(i+(dim*dim)-1)), (dim,dim))
#     elseif rank == 3
#         return reshape(view(input, i:(i+(dim*dim*dim)-1)), (dim, dim,dim))
#     else
#         @warn "tensor_view for rank > 3 is a general implementation that needs allocations!"
#         return reshape(view(input, i:(i+(dim^rank)-1)),ntuple(i->dim,rank))
#     end
# end


"""
````
function tmul!(y,A,x,α=1.0,β=0.0)
````

Combined inplace  matrix-vector multiply-add ``A^T x α + y β``.
The result is stored in `y` by overwriting it.  Note that `y` must not be
aliased with either `A` or `x`.

"""
function tmul!(y, A, x, α = 1.0, β = 0.0)
    for i in eachindex(y)
        y[i] *= β
        for j in eachindex(x)
            y[i] += α * A[j, i] * x[j]
        end
    end
    return
end

"""
````
function tmul!(y::AbstractVector{T}, A::AbstractMatrix{T}, x::AbstractVector{T}, α=1.0, β=0.0) where {T<:AbstractFloat}
````

Overload of the generic function for types supported by 
`LinearAlgebra.BLAS.gemv!` to avoid slow run times for large inputs.
"""
function tmul!(
        y::AbstractVector{T},
        A::AbstractMatrix{T},
        x::AbstractVector{T},
        α = 1.0,
        β = 0.0
    ) where {T <: AbstractFloat}
    return LinearAlgebra.BLAS.gemv!('T', α, A, x, β, y)
end
