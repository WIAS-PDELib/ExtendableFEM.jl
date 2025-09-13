struct MeanValueRestriction{T, Tv} <: AbstractRestriction
    unknown::Unknown
    value::T
    linear_operator
    parameters::Dict{Symbol, Any}
end


"""
$(TYPEDSIGNATURES)

Construct a mean value restriction for a scalar unknown.

This restriction enforces that the mean value of the given unknown
(or rather the testing of an underlying LinearOperator) matches a specified value.  
Optionally, an operator can be applied to the unknown before computing the mean,
and a custom kernel function for the underlying `LinearOperator` can be provided.

# Arguments
- `u::Unknown`: The unknown whose mean value is to be restricted.
- `kernel`: Kernel function for the linear operator (default: `ExtendableFEMBase.constant_one_kernel`).
- `value::T`: The target mean value (default: `0`).
- `operator`: Operator to apply to the unknown before restriction (default: `Identity`).
- `Tv`: Value type for the restriction (default: `Float64`).

"""
function MeanValueRestriction(u::Unknown; kernel = ExtendableFEMBase.constant_one_kernel, value::T = 0, operator = Identity, Tv = Float64) where {T}
    linear_operator = LinearOperator(kernel, [(u, operator)]; store = true, Tv)
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
