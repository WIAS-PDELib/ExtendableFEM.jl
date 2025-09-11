struct MeanValueRestriction{T, Tv} <: AbstractRestriction
    unknown::Unknown
    value::T
    linear_operator
    parameters::Dict{Symbol, Any}
end


"""
$(TYPEDSIGNATURES)

Restrict the mean value of a scalar unknown to a certain value.

An additional operator can be applied to the unknown before the mean value is restricted.
"""
function MeanValueRestriction(u::Unknown; value::T = 0, operator = Identity, Tv = Float64) where {T}
    linear_operator = LinearOperator([(u, operator)]; store = true, Tv)
    return MeanValueRestriction{T, Tv}(u, value, linear_operator, Dict{Symbol, Any}(:name => "MeanValueRestriction"))
end


function assemble!(R::MeanValueRestriction{T, Tv}, sol, SC; kwargs...) where {T, Tv}

    n = length(SC.b.entries)

    # type and shape of the rhs
    b = deepcopy(SC.b)
    fill!(b.entries, 0)

    assemble!(nothing, b, sol, R.linear_operator, SC; assemble_rhs = true, kwargs...)

    R.parameters[:matrix] = reshape(b.entries, n, 1)
    R.parameters[:rhs] = Tv[R.value]

    return nothing
end
