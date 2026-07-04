function structure_constants(ip::AbstractInnerProduct, B::MatrixBasis)
    n = length(B)
    T = promote_type(eltype(B[1]), ComplexF64)
    f = Array{T,3}(undef, n, n, n)

    for i in 1:n
        for j in 1:n
            f[:, i, j] = bracket_coefficients(ip, B, i, j)
        end
    end
    return f
end

function structure_constants(B::MatrixBasis)
    return structure_constants(DEFAULT_INNER_PRODUCT, B)
end
