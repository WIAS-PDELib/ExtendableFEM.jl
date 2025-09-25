"""
	AbstractRestriction

Root type for all restrictions
"""
abstract type AbstractRestriction end


function Base.show(io::IO, R::AbstractRestriction)
    print(io, "AbstractRestriction")
    return nothing
end

# informs solver when operator needs reassembly in a time dependent setting
function is_timedependent(R::AbstractRestriction)
    return false
end

function assemble!(R::AbstractRestriction, SC; kwargs...)
    ## assembles internal restriction matrix in R
    return nothing
end

restriction_matrix(R::AbstractRestriction) = R.parameters[:matrix]
restriction_rhs(R::AbstractRestriction) = R.parameters[:rhs]
fixed_dofs(R::AbstractRestriction) = R.parameters[:fixed_dofs]
