struct LinearFunctionalRestriction{T, Tv} <: AbstractRestriction
    value::T
    linear_operator
    parameters::Dict{Symbol, Any}
end


"""
$(TYPEDSIGNATURES)

Construct a linear functional restriction for a finite element unknown from a given `LinearOperator`.

# Arguments
- `O::LinearOperator{Tv}`: The linear operator representing the functional to be restricted.
- `value::T`: The target value for the linear functional (default: `0`).

"""
function LinearFunctionalRestriction(O::LinearOperator{Tv}; value::T = 0) where {Tv, T}
    return LinearFunctionalRestriction{T, Tv}(value, O, Dict{Symbol, Any}(:name => "LinearFunctionalRestriction"))
end


"""
$(TYPEDSIGNATURES)

Construct a zero mean value restriction for a finite element unknown.

This restriction enforces that the mean value of the given unknown (or of an operator applied to it) vanishes.  
Internally, it creates a `LinearOperator` using the provided kernel and operator, and wraps it in a `LinearFunctionalRestriction`.

# Arguments
- `u::Unknown`: The unknown whose mean value is to be restricted.
- `kernel`: Kernel function for the linear operator (default: `ExtendableFEMBase.constant_one_kernel`).
- `operator`: Operator to apply to the unknown before restriction (default: `Identity`).
- `Tv`: Value type for the restriction (default: `Float64`).

"""
function ZeroMeanValueRestriction(u::Unknown; kernel = ExtendableFEMBase.constant_one_kernel, operator = Identity, Tv = Float64)
    linear_operator = LinearOperator(kernel, [(u, operator)]; store = true, Tv)
    return LinearFunctionalRestriction{Int64, Tv}(0, linear_operator, Dict{Symbol, Any}(:name => "MeanValueRestriction"))
end


"""
$(TYPEDSIGNATURES)

Construct a mass/integral restriction for a finite element unknown.

This restriction enforces that the integral of the given unknown (or of an operator applied to it) attains a given value.  
Internally, it creates a `LinearOperator` using the provided kernel and operator, and wraps it in a `LinearFunctionalRestriction`.

# Arguments
- `u::Unknown`: The unknown whose mean value is to be restricted.
- `kernel`: Kernel function for the linear operator (default: `ExtendableFEMBase.constant_one_kernel`).
- `value::T`: The target integral/mass value (default: `0`).
- `operator`: Operator to apply to the unknown before restriction (default: `Identity`).
- `Tv`: Value type for the restriction (default: `Float64`).

"""
function MassRestriction(u::Unknown; kernel = ExtendableFEMBase.constant_one_kernel, value::T = 0, operator = Identity, Tv = Float64) where {T}
    linear_operator = LinearOperator(kernel, [(u, operator)]; store = true, Tv)
    return LinearFunctionalRestriction{T, Tv}(value, linear_operator, Dict{Symbol, Any}(:name => "MeanValueRestriction"))
end


function assemble!(R::LinearFunctionalRestriction{T, Tv}, sol, SC; kwargs...) where {T, Tv}

    n = length(SC.b.entries)

    b = deepcopy(SC.b)
    fill!(b.entries, 0.0)

    assemble!(nothing, b, sol, R.linear_operator, SC; assemble_rhs = true, kwargs...)

    R.parameters[:matrix] = reshape(b.entries, n, 1)
    R.parameters[:rhs] = Tv[R.value]
    R.parameters[:fixed_dofs] = []

    return nothing
end
