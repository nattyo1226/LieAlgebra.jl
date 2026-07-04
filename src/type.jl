function check_empty(elements::AbstractVector{<:AbstractMatrix})
    isempty(elements) && throw(ArgumentError("Elements cannot be empty."))
    return true
end

function check_same_size(elements::AbstractVector{<:AbstractMatrix})
    s = size(first(elements))
    for e in elements
        size(e) == s || throw(DimensionMismatch("All matrices must have the same size."))
    end
    return true
end

function check_square(elements::AbstractVector{<:AbstractMatrix})
    for e in elements
        size(e, 1) == size(e, 2) || throw(DimensionMismatch("All matrices must be square."))
    end
    return true
end

struct MatrixBasis{M<:AbstractMatrix}
    elements::Vector{M}

    function MatrixBasis(elements::AbstractVector{M}) where {M<:AbstractMatrix}
        check_empty(elements)
        check_same_size(elements)
        check_square(elements)

        new{M}(collect(elements))
    end
end
