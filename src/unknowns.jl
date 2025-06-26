"""
$(TYPEDEF)

A structure representing a problem "unknown" (e.g., a field variable) with associated metadata.

# Fields
- `name::String`: name of the unknown (used in printouts and messages).
- `identifier::IT`: identifier for operator assignment and indexing.
- `parameters::Dict{Symbol, Any}`: Dictionary of additional properties (e.g., dimension, symbols for ansatz/test functions, ...).

# Type Parameters
- `IT`: Type of the identifier (commonly `Symbol`).

# Usage
Create an `Unknown` using the provided constructor, optionally specifying the name, identifier, and additional keyword arguments for parameters.

# Example
```julia
u = Unknown("velocity"; dimension=2, symbol_ansatz="u", symbol_test="v")
```
"""
mutable struct Unknown{IT}
    """
    The name of the unknown used for printout messages.
    """
    name::String
    """
    The identifier of the unknown used for assignments to operators.
    """
    identifier::IT
    """
    Further properties of the unknown can be stored in a Dict, see constructor.
    """
    parameters::Dict{Symbol, Any}
end


default_unknown_kwargs() = Dict{Symbol, Tuple{Any, String}}(
    :dimension => (nothing, "dimension of the unknown"),
    :symbol_ansatz => (nothing, "symbol for ansatz functions of this unknown in printouts"),
    :symbol_test => (nothing, "symbol for test functions of this unknown in printouts"),
    :algebraic_constraint => (nothing, "is this unknown an algebraic constraint?"),
)

test_function(u::Unknown) = u.parameters[:symbol_test] === nothing ? "X($(u.identifier))" : u.parameters[:symbol_test]
ansatz_function(u::Unknown) = u.parameters[:symbol_ansatz] === nothing ? u.identifier : u.parameters[:symbol_ansatz]


"""
$(TYPEDSIGNATURES)

Construct an `Unknown` representing a finite element field variable or problem unknown.

# Arguments
- `u::String`: Name of the unknown.
- `identifier`: (optional) Symbolic or integer identifier for operator assignment and indexing (default: `Symbol(u)`).
- `name`: (optional) name in printouts (default: `u`).
- `kwargs...`: Additional keyword arguments to set or override entries in the `parameters` dictionary (see below).

# Keyword Arguments
$(_myprint(default_unknown_kwargs()))

# Example
```julia
u = Unknown("velocity"; dimension=2, symbol_ansatz="u", symbol_test="v")

"""
function Unknown(u::String; identifier = Symbol(u), name = u, kwargs...)
    parameters = Dict{Symbol, Any}(k => v[1] for (k, v) in default_unknown_kwargs())
    _update_params!(parameters, kwargs)
    return Unknown(name, identifier, parameters)
end

function Base.show(io::IO, u::Unknown)
    return print(io, "$(u.identifier) ($(u.name))")
end

## remapping of all function operators
FO(u) = (u, FO)

"""
	jump(o:Tuple{Union{Unknown, Int}, StandardFunctionOperator})

alias for (o[1], Jump{o[2]})
"""
jump(o::Tuple{Union{Unknown, Int}, DataType}) = (o[1], Jump{o[2]})

"""
	average(o:Tuple{Union{Unknown, Int}, StandardFunctionOperator})

alias for (o[1], Average{o[2]})
"""
average(o::Tuple{Union{Unknown, Int}, DataType}) = (o[1], Average{o[2]})

"""
	this(o:Tuple{Union{Unknown, Int}, StandardFunctionOperator})

alias for (o[1], Left{o[2]})
"""
this(o::Tuple{Union{Unknown, Int}, DataType}) = (o[1], Left{o[2]})

"""
	other(o:Tuple{Union{Unknown, Int}, StandardFunctionOperator})

alias for (o[1], Right{o[2]})
"""
other(o::Tuple{Union{Unknown, Int}, DataType}) = (o[1], Right{o[2]})

## some aliases

"""
	id(u)

alias for (u, Identity)
"""
id(u) = (u, Identity)

"""
	curl1(u)

alias for (u, CurlScalar)
"""
curl1(u) = (u, CurlScalar)

"""
	curl2(u)

alias for (u, Curl2D)
"""
curl2(u) = (u, Curl2D)


"""
	curl3(u)

alias for (u, Curl3D)
"""
curl3(u) = (u, Curl3D)

"""
	grad(u)

alias for (u, Gradient)
"""
grad(u) = (u, Gradient)

"""
	laplace(u)

alias for (u, Laplacian)
"""
laplace(u) = (u, Laplacian)

"""
	Δ(u)

alias for (u, Laplacian)
"""
Δ(u) = (u, Laplacian)

"""
	hessian(u)

alias for (u, Hessian)
"""
hessian(u) = (u, Hessian)

"""
	div(u)

alias for (u, Divergence)
"""
Base.div(u) = (u, Divergence)

"""
	normalflux(u)

alias for (u, NormalFlux)
"""
normalflux(u) = (u, NormalFlux)

"""
	tangentialflux(u)

alias for (u, TangentFlux)
"""
tangentialflux(u) = (u, TangentFlux)

"""
	tangentialgrad(u)

alias for (u, TangentialGradient)
"""
tangentialgrad(u) = (u, TangentialGradient)

"""
	symgrad_voigt(u, factor)

alias for (u, SymmetricGradient{factor}) in Voigt notation.

The `factor` is a real number applied to the (summed) off-diagonal entries.
In the current implementation in ExtendableFEMBase, factor = 0.5 represents the Voigt strain mapping [σ₁₁ σ₂₂ σ₃₃ σ₁₃ σ₂₃ σ₁₂],
while factor = 1.0 represents the Voigt strain mapping [ε₁₁ ε₂₂ ε₃₃ 2ε₁₃ 2ε₂₃ 2ε₁₂].
"""
symgrad_voigt(u, factor) = (u, SymmetricGradient{factor})

"""
	εV(u, factor)

unicode alias for [`symgrad_voigt(u, factor)`](@ref).
"""
εV(u, factor) = symgrad_voigt(u, factor)

"""
	grid(u)

alias for (u, "grid") (triggers gridplot in plot)
"""
grid(u) = (u, "grid")

"""
	dofgrid(u)

alias for (u, "dofgrid") (triggers gridplot for the dofgrid in plot)
"""
dofgrid(u) = (u, "dofgrid")


"""
	apply(u, FO::Type{<:AbstractFunctionOperator})

alias for (u, FO)
"""
apply(u, FO::Type{<:AbstractFunctionOperator}) = (u, FO)
