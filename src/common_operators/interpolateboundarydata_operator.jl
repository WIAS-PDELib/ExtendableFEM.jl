mutable struct InterpolateBoundaryData{UT, DFT} <: AbstractOperator
    u::UT
    data::DFT
    bdofs::Array{Int, 1}
    bfaces::Array{Int, 1}
    FES::Any
    bddata::Any
    assembler::Any
    parameters::Dict{Symbol, Any}
end

"""
````
fixed_dofs(O::InterpolateBoundaryData)
````

returns the fixed degrees of freedoms of O
"""
fixed_dofs(O::InterpolateBoundaryData) = O.bdofs

"""
````
fixed_vals(O::InterpolateBoundaryData)
````

returns the currently assembled values for the fixed degrees of freedom of O
"""
fixed_vals(O::InterpolateBoundaryData) = view(O.bddata.entries, O.bdofs)

default_bndop_kwargs() = Dict{Symbol, Tuple{Any, String}}(
    :penalty => (1.0e30, "penalty for fixed degrees of freedom"),
    :name => ("BoundaryData", "name for operator used in printouts"),
    :bonus_quadorder => (0, "additional quadrature order added to the quadorder chosen by the interpolator"),
    :params => (nothing, "array of parameters that should be made available in qpinfo argument of kernel function"),
    :regions => ([], "subset of regions where operator should be assembly only"),
    :plot => (false, "plot unicode plot of boundary data into terminal when assembled"),
    :verbosity => (0, "verbosity level"),
)

# informs solver in which blocks the operator assembles to
function dependencies_when_linearized(O::InterpolateBoundaryData)
    return O.u
end

# informs solver when operator needs reassembly in a time dependent setting
function is_timedependent(O::InterpolateBoundaryData)
    return O.parameters[:time_dependent]
end

function Base.show(io::IO, O::InterpolateBoundaryData)
    dependencies = dependencies_when_linearized(O)
    print(io, "$(O.parameters[:name])($(ansatz_function(dependencies)))")
    return nothing
end

"""
````
function InterpolateBoundaryData(u, data!::Function; kwargs...)
````

When assembled, the unknown u of the Problem will be penalized
to match the standard interpolation of the provided data! function.
The header of this function needs to be conform to the interface

	data!(result, qpinfo)

where qpinfo allows to access information at the current quadrature point,
e.g. qpinfo.x provides the global coordinates of the quadrature/evaluation point.

Keyword arguments:
$(_myprint(default_bndop_kwargs()))

"""
function InterpolateBoundaryData(u, data = nothing; kwargs...)
    parameters = Dict{Symbol, Any}(k => v[1] for (k, v) in default_bndop_kwargs())
    _update_params!(parameters, kwargs)
    return InterpolateBoundaryData{typeof(u), typeof(data)}(u, data, zeros(Int, 0), zeros(Int, 0), nothing, nothing, nothing, parameters)
end


"""
````
assemble!(A, b, sol, O::InterpolateBoundaryData{UT}, SC::SolverConfiguration; kwargs...)
````

assembles the correct boundary values for O by interpolating the boundary data with
the current kwargs (where e.g. time and params might have changed).
"""
function assemble!(A, b, sol, O::InterpolateBoundaryData{UT}, SC::SolverConfiguration; kwargs...) where {UT}
    if UT <: Integer
        ind = O.u
        ind_sol = ind
    elseif UT <: Unknown
        ind = get_unknown_id(SC, O.u)
        ind_sol = findfirst(==(O.u), sol.tags)
    end
    return assemble!(O, b[ind].FES; ind = ind, offset = SC.offsets[ind], kwargs...)
end

function assemble!(O::InterpolateBoundaryData, FES = O.FES; time = 0, offset = 0, kwargs...)
    regions = O.parameters[:regions]
    bdofs::Array{Int, 1} = O.bdofs
    bfaces::Array{Int, 1} = O.bfaces
    if O.FES !== FES
        bddata = FEVector(FES)
        xgrid = FES.dofgrid
        Ti = eltype(xgrid[CellNodes])
        bfacedofs::Adjacency{Ti} = FES[BFaceDofs]
        bfacefaces = xgrid[BFaceFaces]
        bfaceregions = xgrid[BFaceRegions]
        nbfaces = num_sources(bfacedofs)
        ndofs4bface = max_num_targets_per_source(bfacedofs)
        bdofs = []
        bfaces = []
        for bface in 1:nbfaces
            if bfaceregions[bface] in regions
                for k in 1:ndofs4bface
                    dof = bfacedofs[k, bface] + offset
                    push!(bdofs, dof)
                end
                push!(bfaces, bfacefaces[bface])
            end
        end
        unique!(bdofs)
        O.bdofs = bdofs
        O.bddata = bddata
        O.bfaces = bfaces
        if O.parameters[:verbosity] > 0
            @info "$(O.parameters[:name]) : penalizing $(length(O.bdofs)) dofs of '$(O.u.name)' (ON_BFACES, regions = $(O.parameters[:regions]))"
        end
    end
    time = @elapsed begin
        bddata = O.bddata
        data = O.data
        bfaces = O.bfaces
        if FES.broken
            FEType = eltype(FES)
            xgrid = FES.dofgrid
            FESc = FESpace{FEType}(xgrid)
            Targetc = FEVector(FESc)
            interpolate!(Targetc[1], FESc, ON_FACES, data; items = bfaces, time = time, params = O.parameters[:params], bonus_quadorder = O.parameters[:bonus_quadorder])
            bfacedofs = FES[BFaceDofs]
            bfacedofs_c = FESc[BFaceDofs]
            dof::Int = 0
            dofc::Int = 0
            for bface in 1:nbfaces
                for k in 1:num_targets(bfacedofs, bface)
                    dof = bfacedofs[k, bface]
                    dofc = bfacedofs_c[k, bface]
                    bddata.entries[dof] = Targetc.entries[dofc]
                end
            end
        else
            interpolate!(bddata[1], ON_BFACES, data; time = time, params = O.parameters[:params], bonus_quadorder = O.parameters[:bonus_quadorder])
        end
        if O.parameters[:plot]
            println(stdout, unicode_scalarplot(bddata[1]; title = "boundary data for $(O.u)"))
        end
    end
    return if O.parameters[:verbosity] > 0
        @info "$(O.parameters[:name]) : assembly took $time s"
    end
end


"""
````
apply!(U::FEVectorBlock, O::InterpolateBoundaryData; offset = 0, kwargs...)
````

applies the boundary data of O to U, i.e., sets the boundary dofs to the correct values
that have been computed during the last assemble! call.
"""
function apply!(U::FEVectorBlock, O::InterpolateBoundaryData; offset = 0, kwargs...)
    bddata = O.bddata
    bdofs = O.bdofs
    for dof in bdofs
        U[dof - offset] = bddata.entries[dof - offset]
    end
    return
end

"""
````
apply_penalties!(A, b, sol, O::InterpolateBoundaryData{UT}, SC::SolverConfiguration; kwargs...)
````

modifies the linear system A|b such that the boundary dofs are penalized and attain the
correct values from the last assemble! call of O. Also applies the correct values to sol.
"""
function apply_penalties!(A, b, sol, O::InterpolateBoundaryData{UT}, SC::SolverConfiguration; assemble_matrix = true, assemble_rhs = true, assemble_sol = true, kwargs...) where {UT}
    time = @elapsed begin
        if UT <: Integer
            ind = O.u
            ind_sol = ind
        elseif UT <: Unknown
            ind = get_unknown_id(SC, O.u)
            ind_sol = findfirst(==(O.u), sol.tags)
        end
        offset = SC.offsets[ind]
        bddata = O.bddata
        bdofs = O.bdofs
        penalty = O.parameters[:penalty]
        if assemble_matrix
            AE = A.entries
            for dof in bdofs
                rawupdateindex!(AE, (a, b) -> b, penalty, dof, dof, 1)
            end
            flush!(AE)
        end
        if assemble_rhs
            BE = b.entries
            for dof in bdofs
                BE[dof] = penalty * bddata.entries[dof - offset]
            end
        end
        if assemble_sol
            apply!(sol[ind_sol], O; offset = offset)
        end
    end
    return if O.parameters[:verbosity] > 1
        @info "$(O.parameters[:name]) : applying penalties took $time s"
    end
end
