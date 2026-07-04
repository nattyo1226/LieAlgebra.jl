commutator(A::AbstractMatrix, B::AbstractMatrix) = A * B - B * A

function bracket(B::MatrixBasis, i::Integer, j::Integer)
    return commutator(B[i], B[j])
end

function bracket_coefficients(ip::AbstractInnerProduct, B::MatrixBasis, i::Integer, j::Integer)
    return coefficients(ip, B, bracket(B, i, j))
end

function bracket_coefficients(B::MatrixBasis, i::Integer, j::Integer)
    return bracket_coefficients(DEFAULT_INNER_PRODUCT, B, i, j)
end
