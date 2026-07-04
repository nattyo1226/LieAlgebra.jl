function normalize_matrix(ip::AbstractInnerProduct, A::AbstractMatrix; atol::Real=1e-10)
    nrm2 = real(ip(A, A))
    nrm2 > atol || throw(ArgumentError("Cannot normalize a zero-norm matrix."))
    return A / sqrt(nrm2)
end

normalize_matrix(A::AbstractMatrix; atol::Real=1e-10) = normalize_matrix(DEFAULT_INNER_PRODUCT, A; atol=atol)

function orthogonalize(ip::AbstractInnerProduct, B::MatrixBasis; atol::Real=1e-10)
    orthogonal = Vector{Matrix{ComplexF64}}()

    for A in B
        V = Matrix{ComplexF64}(A)

        for Q in orthogonal
            V -= ip(Q, V) * Q
        end

        if real(ip(V, V)) > atol
            push!(orthogonal, normalize_matrix(ip, V; atol=atol))
        end
    end

    length(orthogonal) == length(B) || throw(ArgumentError("Basis contains linearly dependent elements."))

    return MatrixBasis(orthogonal)
end

orthogonalize(B::MatrixBasis; atol::Real=1e-10) = orthogonalize(DEFAULT_INNER_PRODUCT, B; atol=atol)

function orthogonalize(ip::AbstractInnerProduct, D::Decomposition; atol::Real=1e-10)
    ideals = [
        orthogonalize(ip, B; atol=atol)
        for B in D.ideals
    ]

    return Decomposition(ideals)
end

orthogonalize(D::Decomposition; atol::Real=1e-10) = orthogonalize(DEFAULT_INNER_PRODUCT, D; atol=atol)

function is_orthonormal(ip::AbstractInnerProduct, B::MatrixBasis; atol::Real=1e-10)
    G = gram_matrix(ip, B)
    return norm(G - I(length(B))) < atol
end

is_orthonormal(B::MatrixBasis; atol::Real=1e-10) = is_orthonormal(DEFAULT_INNER_PRODUCT, B; atol=atol)
