function _check_empty(As::AbstractVector{<:AbstractMatrix})
    isempty(As) && throw(ArgumentError("Input vector must not be empty."))
    return true
end

function _check_same_square(As::AbstractVector{<:AbstractMatrix})
    n = size(first(As), 1)
    for A in As
        size(A) == (n, n) || throw(DimensionMismatch("All matrices must have the same size."))
    end
    return true
end

struct MatrixBasis{M<:AbstractMatrix}
    elements::Vector{M}

    function MatrixBasis(elements::AbstractVector{M}) where {M<:AbstractMatrix}
        _check_empty(elements)
        _check_same_square(elements)

        new{M}(collect(elements))
    end
end

struct Decomposition
    ideals::Vector{MatrixBasis}
end
