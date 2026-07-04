function killing_form(ip::AbstractInnerProduct, B::MatrixBasis)
    ads = adjoint_matrices(ip, B)

    n = length(B)
    K = Matrix{ComplexF64}(undef, n, n)

    for i in 1:n
        for j in 1:n
            K[i, j] = ip(ads[i], ads[j])
        end
    end

    return K
end

function killing_form(B::MatrixBasis)
    return killing_form(DEFAULT_INNER_PRODUCT, B)
end
