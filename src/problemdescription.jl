"""
$(TYPEDEF)

A structure representing a finite element problem description, including unknowns and operators.

# Fields
- `name::String`: name of the problem (used in printouts and logs). Default: `"My Problem"`.
- `unknowns::Vector{Unknown}`: List of unknowns involved in the problem.
- `operators::Vector{AbstractOperator}`: List of operators (e.g., bilinear forms, linear forms, boundary conditions) that define the problem.

# Usage
Create a `ProblemDescription` using the default constructor, then assign unknowns and operators using `assign_unknown!` and `assign_operator!`.

"""
mutable struct ProblemDescription
    """
    The name of the problem used for printout messages. Default: "My Problem"
    """
    name::String
    """
    A vector of Unknowns that are involved in the problem.
    """
    unknowns::Array{Unknown, 1}
    """
    A vector of operators that are involved in the problem.
    """
    operators::Array{AbstractOperator, 1}
    #reduction_operators::Array{AbstractReductionOperator,1}

    """
    A vector of Lagrange restrictions that are involved in the problem.
    """
    restrictions::Vector{AbstractRestriction}
end

"""
$(TYPEDSIGNATURES)

Create an empty `ProblemDescription` with the given name.

# Arguments
- `name::String`: (optional) Name of the problem (used in printouts and logs). Default: `"My problem"`.

# Returns
- A `ProblemDescription` instance with the specified name and empty unknown/operator lists.

"""
function ProblemDescription(name = "My problem")
    return ProblemDescription(name, [], [], [])
end


"""
$(TYPEDSIGNATURES)

Add an `Unknown` to a `ProblemDescription`, if not already present.

# Arguments
- `PD::ProblemDescription`: The problem description to which the unknown should be added.
- `u::Unknown`: The unknown to add.

# Returns
- The index (position) of the unknown in the `unknowns` array of `PD`.

# Notes
- If the unknown is already present, a warning is issued and its existing position is returned.


"""
function assign_unknown!(PD::ProblemDescription, u::Unknown)
    if u in PD.unknowns
        @warn "This unknown was already assigned to the problem description! Ignoring this call."
        return find(==(u), PD.unknowns)
    else
        push!(PD.unknowns, u)
        return length(PD.unknowns)
    end
end


"""
$(TYPEDSIGNATURES)

Adds an operator to a `ProblemDescription`.

# Arguments
- `PD::ProblemDescription`: The problem description to which the operator should be added.
- `o::AbstractOperator`: The operator to add.

# Returns
- The index (position) of the operator in the `operators` array of `PD`.

"""
function assign_operator!(PD::ProblemDescription, o::AbstractOperator)
    push!(PD.operators, o)
    return length(PD.operators)
end

"""
$(TYPEDSIGNATURES)

Adds a restrction to a `ProblemDescription`.

# Arguments
- `PD::ProblemDescription`: The problem description to which the operator should be added.
- `r::AbstractRestriction`: The restriction to add.

# Returns
- The index (position) of the restriction in the `restrictions` array of `PD`.
"""
function assign_restriction!(PD::ProblemDescription, r::AbstractRestriction)
    push!(PD.restrictions, r)
    return length(PD.restrictions)
end


"""
$(TYPEDSIGNATURES)

Replace an operator in a `ProblemDescription` at a specified position.

# Arguments
- `PD::ProblemDescription`: The problem description in which to replace the operator.
- `j::Int`: The index (position) of the operator to replace.
- `o::AbstractOperator`: The new operator to insert.

"""
function replace_operator!(PD::ProblemDescription, j, o::AbstractOperator)
    PD.operators[j] = o
    return nothing
end

#function assign_reduction!(PD::ProblemDescription, u::AbstractReductionOperator)
#    push!(PD.reduction_operators, u)
#end

function Base.show(io::IO, PD::ProblemDescription)
    println(io, "\nPDE PROBLEM DESCRIPTION")
    println(io, "  Name: $(PD.name)")
    println(io, "  Unknowns ($(length(PD.unknowns))):")
    if isempty(PD.unknowns)
        println(io, "    (none)")
    else
        for (i, u) in enumerate(PD.unknowns)
            println(io, "    [$i] $u")
        end
    end
    println(io, "  Operators ($(length(PD.operators))):")
    if isempty(PD.operators)
        println(io, "    (none)")
    else
        for (i, o) in enumerate(PD.operators)
            println(io, "    [$i] $o")
        end
    end
    return
end
