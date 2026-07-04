function adjoint_matrix(ip::AbstractInnerProduct, B::MatrixBasis, i::Integer)
    n = length(B)

    ad = Matrix{ComplexF64}(undef, n, n)

    for j in 1:n
        ad[:, j] .= bracket_coefficients(ip, B, i, j)
    end

    return ad
end

function adjoint_matrix(B::MatrixBasis, i::Integer)
    return adjoint_matrix(DEFAULT_INNER_PRODUCT, B, i)
end

function adjoint_matrices(ip::AbstractInnerProduct, B::MatrixBasis)
    return [adjoint_matrix(ip, B, i) for i in eachindex(B)]
end

function adjoint_matrices(B::MatrixBasis)
    return adjoint_matrices(DEFAULT_INNER_PRODUCT, B)
end
