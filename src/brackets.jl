commutator(A::AbstractMatrix, B::AbstractMatrix) = A * B - B * A

function lie_bracket(A::AbstractMatrix, B::AbstractMatrix)
    return -im * commutator(A, B)
end

function bracket(B::MatrixBasis, i::Integer, j::Integer)
    return lie_bracket(B[i], B[j])
end

function bracket_coefficients(ip::AbstractInnerProduct, B::MatrixBasis, i::Integer, j::Integer)
    return coefficients(ip, B, bracket(B, i, j))
end

function bracket_coefficients(B::MatrixBasis, i::Integer, j::Integer)
    return bracket_coefficients(DEFAULT_INNER_PRODUCT, B, i, j)
end
