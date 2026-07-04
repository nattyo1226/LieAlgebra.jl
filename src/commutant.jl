function commutant_constraint(A::AbstractMatrix)
    n = size(A, 1)
    size(A, 2) == n || throw(DimensionMismatch("Matrix must be square."))

    In = Matrix{eltype(A)}(I, n, n)
    return kron(transpose(A), In) - kron(In, A)
end

function commutant_matrix(As::AbstractVector{<:AbstractMatrix})
    _check_empty(As)
    _check_same_square(As)

    return vcat(commutant_constraint.(As)...)
end

function commutant(As::AbstractVector{<:AbstractMatrix}; atol::Real=1e-10)
    M = commutant_matrix(As)
    ns = nullspace(M; atol=atol)

    n = size(first(As), 1)

    return [reshape(ns[:, i], (n, n)) for i in axes(ns, 2)]
end

function commutant(ip::AbstractInnerProduct, B::MatrixBasis; atol::Real=1e-10)
    As = adjoint_matrices(ip, B)
    return commutant(As; atol=atol)
end

function commutant(B::MatrixBasis; atol::Real=1e-10)
    return commutant(DEFAULT_INNER_PRODUCT, B; atol=atol)
end
