struct BoundaryDataRestriction{BOT} <: AbstractRestriction
    unknown::Unknown
    boundary_operator::BOT
    parameters::Dict{Symbol, Any}
end


"""
$(TYPEDSIGNATURES)

Construct a boundary data restriction for an unknown.

This restriction enforces prescribed boundary values for the given unknown, either as homogeneous data (similar to a HomogeneousBoundaryData operator)
or interpolated from a provided data function (similar to an InterpolateBoundaryOperator).
    
# Arguments
- `u::Unknown`: The unknown whose boundary values are to be restricted.
- `data`: A function with signature `data!(result, qpinfo)` that computes the desired boundary values at each quadrature point. The `qpinfo` argument provides information such as global coordinates (`qpinfo.x`). If no
data is given, homogeneous boundary data is assembled
- `kwargs...`: Additional keyword arguments passed to the boundary operator constructors.

"""
function BoundaryDataRestriction(u::Unknown, data = nothing; kwargs...)
    if data === nothing
        boundary_operator = HomogeneousBoundaryData(u; kwargs...)
    else
        boundary_operator = InterpolateBoundaryData(u, data; kwargs...)
    end
    return BoundaryDataRestriction(u, boundary_operator, Dict{Symbol, Any}(:name => "BoundaryDataRestriction"))
end


function assemble!(R::BoundaryDataRestriction, sol, SC; kwargs...)

    ## assemble boundary operator (-> computes fixed dofs and vals)
    assemble!(nothing, SC.b, sol, R.boundary_operator, SC; kwargs...)
    fixeddofs = fixed_dofs(R.boundary_operator)
    fixedvals = fixed_vals(R.boundary_operator)
    nvals = length(fixeddofs)

    ## assign fixed dofs and vals to restriction
    n = length(SC.b.entries)
    R.parameters[:matrix] = sparse(fixeddofs, 1:nvals, ones(nvals), n, nvals)
    R.parameters[:rhs] = fixedvals

    return nothing
end
