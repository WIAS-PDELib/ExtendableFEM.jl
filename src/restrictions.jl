"""
	AbstractRestriction

Root type for all restrictions
"""
abstract type AbstractRestriction end


function Base.show(io::IO, R::AbstractRestriction)
    if haskey(R.parameters, :name)
        print(io, "$(R.parameters[:name])")
        if haskey(R.parameters, :name)
            print(io, " for $(ansatz_function(R.parameters[:unknown]))")
        end
    else
        print(io, "Restriction of type $(typeof(R))")
    end
    return nothing
end

# informs solver when operator needs reassembly in a time dependent setting
function is_timedependent(R::AbstractRestriction)
    return false
end

function assemble!(R::AbstractRestriction, SC; kwargs...)
    ## assembles internal restriction matrix in R
    @warn "assemble! not implemented for $(typeof(R))"
    return nothing
end

restriction_matrix(R::AbstractRestriction) = R.parameters[:matrix]
restriction_rhs(R::AbstractRestriction) = R.parameters[:rhs]
