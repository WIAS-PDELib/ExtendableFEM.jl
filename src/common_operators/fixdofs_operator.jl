##################
### FIXED DOFS ###
##################

mutable struct FixDofs{UT, AT, VT} <: AbstractOperator
    u::UT
    dofs::AT
    offset::Int
    vals::VT
    assembler::Any
    parameters::Dict{Symbol, Any}
end

fixed_dofs(O::FixDofs) = O.dofs .+ O.offset
fixed_vals(O::FixDofs) = O.vals

default_fixdofs_kwargs() = Dict{Symbol, Tuple{Any, String}}(
    :penalty => (1.0e30, "penalty for fixed degrees of freedom"),
    :name => ("FixDofs", "name for operator used in printouts"),
    :verbosity => (0, "verbosity level"),
)

# informs solver in which blocks the operator assembles to
function dependencies_when_linearized(O::FixDofs)
    return O.u
end

function Base.show(io::IO, O::FixDofs)
    dependencies = dependencies_when_linearized(O)
    print(io, "$(O.parameters[:name])($(ansatz_function(dependencies)), ndofs = $(length(O.dofs)))")
    return nothing
end

"""
$(TYPEDSIGNATURES)

Construct a `FixDofs` operator to strongly enforce fixed values on specified degrees of freedom (dofs) in a finite element system.

# Arguments
- `u`: The unknown (field variable) or its identifier for which dofs are to be fixed.

# Keyword Arguments
- `dofs`: Vector of dof indices to be fixed (default: empty vector).
- `vals`: Vector of values to assign to the fixed dofs (default: zeros of same length as `dofs`).
$(_myprint(default_fixdofs_kwargs()))

# Returns
- A `FixDofs` operator instance specifying which dofs to fix and their target values.
"""
function FixDofs(u; dofs = [], vals = zeros(Float64, length(dofs)), kwargs...)
    parameters = Dict{Symbol, Any}(k => v[1] for (k, v) in default_fixdofs_kwargs())
    _update_params!(parameters, kwargs)
    @assert length(dofs) == length(vals)
    return FixDofs{typeof(u), typeof(dofs), typeof(vals)}(u, dofs, 0, vals, nothing, parameters)
end

function apply_penalties!(A, b, sol, O::FixDofs{UT}, SC::SolverConfiguration; assemble_matrix = true, assemble_rhs = true, kwargs...) where {UT}
    if UT <: Integer
        ind = O.u
    elseif UT <: Unknown
        ind = get_unknown_id(SC, O.u)
    end
    offset = sol[ind].offset
    dofs = O.dofs
    vals = O.vals
    penalty = O.parameters[:penalty]
    if assemble_matrix
        AE = A.entries
        for j in 1:length(dofs)
            dof = dofs[j] + offset
            rawupdateindex!(AE, (a, b) -> b, penalty, dof, dof, 1)
        end
    end
    if assemble_rhs
        BE = b.entries
        for j in 1:length(dofs)
            dof = dofs[j] + offset
            BE[dof] = penalty * vals[j]
        end
    end
    SE = sol.entries
    for j in 1:length(dofs)
        dof = dofs[j] + offset
        SE[dof] = vals[j]
    end
    return O.offset = offset
end
