abstract type AbstractInnerProduct end

function inner_product(ip::AbstractInnerProduct, ::AbstractMatrix, ::AbstractMatrix)
    throw(ArgumentError("inner_product is not implemented for $(typeof(ip))."))
end

function gram_matrix(ip::AbstractInnerProduct, basis::MatrixBasis)
    n = length(basis)
    G = Matrix{promote_type(eltype(basis[1]), ComplexF64)}(undef, n, n)
    for i in 1:n
        for j in 1:n
            G[i, j] = ip(basis[i], basis[j])
        end
    end
    return G
end

function (ip::AbstractInnerProduct)(A::AbstractMatrix, B::AbstractMatrix)
    size(A) == size(B) || throw(DimensionMismatch("Matrices must have the same size."))
    return inner_product(ip, A, B)
end

struct TraceInnerProduct <: AbstractInnerProduct
    normalize::Bool

    function TraceInnerProduct(; normalize::Bool=true)
        new(normalize)
    end
end

function inner_product(ip::TraceInnerProduct, A::AbstractMatrix, B::AbstractMatrix)
    d = ip.normalize ? size(A, 1) : 1
    return sum(conj.(A) .* B) / d
end

const DEFAULT_INNER_PRODUCT = TraceInnerProduct(normalize=true)
