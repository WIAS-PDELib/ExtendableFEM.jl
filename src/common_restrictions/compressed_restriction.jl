struct CompressedRestriction <: AbstractRestriction
    parameters::Dict{Symbol, Any}
end


"""
    CompressedRestriction(matrix::AbstractMatrix, rhs::AbstractVector)

    Creates a simple linear restriction as a result of the QR compression of multiple restrictions
    in the `solve` call. The compressed matrix and the compressed can be handed over directly in
    the constructor.

    This is 𝑛𝑜𝑡 exported, for internal use only!
"""
function CompressedRestriction(matrix::AbstractMatrix, rhs::AbstractVector)
    return CompressedRestriction(
        Dict{Symbol, Any}(
            :name => "CompressedRestriction",
            :matrix => matrix,
            :rhs => rhs,
            :multiplier => zero(rhs)
        )
    )
end
