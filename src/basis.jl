Base.length(B::MatrixBasis) = length(B.elements)
Base.getindex(B::MatrixBasis, i::Integer) = B.elements[i]
Base.iterate(B::MatrixBasis, state...) = iterate(B.elements, state...)
Base.keys(B::MatrixBasis) = Base.OneTo(length(B))
Base.eachindex(B::MatrixBasis) = Base.OneTo(length(B))

elements(B::MatrixBasis) = B.elements

matrix_size(B::MatrixBasis) = size(first(B.elements))

function basis_matrix(B::MatrixBasis)
    return hcat(vec.(B.elements)...)
end

function LinearAlgebra.rank(B::MatrixBasis; atol::Real=1e-10)
    return LinearAlgebra.rank(basis_matrix(B); atol=atol)
end

function is_independent(B::MatrixBasis; atol::Real=1e-10)
    return rank(B; atol=atol) == length(B)
end

function assert_independent(B::MatrixBasis; atol::Real=1e-10)
    is_independent(B; atol=atol) || throw(ArgumentError("The basis is not linearly independent."))
    return true
end

function coefficients(ip::AbstractInnerProduct, B::MatrixBasis, A::AbstractMatrix)
    size(A) == matrix_size(B) || throw(DimensionMismatch("Matrix size does not match the basis."))

    G = gram_matrix(ip, B)
    b = [ip(e, A) for e in B]

    return G \ b
end

function reconstruct(B::MatrixBasis, coeffs::AbstractVector)
    length(coeffs) == length(B) || throw(DimensionMismatch("Coefficient vector length does not match the basis."))

    return sum(c * b for (c, b) in zip(coeffs, B))
end

function project(ip::AbstractInnerProduct, B::MatrixBasis, A::AbstractMatrix)
    coeffs = coefficients(ip, B, A)
    return reconstruct(B, coeffs)
end

function project(B::MatrixBasis, A::AbstractMatrix)
    return project(DEFAULT_INNER_PRODUCT, B, A)
end

function residual(ip::AbstractInnerProduct, B::MatrixBasis, A::AbstractMatrix)
    return norm(A - project(ip, B, A))
end

function residual(B::MatrixBasis, A::AbstractMatrix)
    return residual(DEFAULT_INNER_PRODUCT, B, A)
end

function in_span(ip::AbstractInnerProduct, B::MatrixBasis, A::AbstractMatrix; atol::Real=1e-10)
    return residual(ip, B, A) < atol
end

function in_span(B::MatrixBasis, A::AbstractMatrix; atol::Real=1e-10)
    return in_span(DEFAULT_INNER_PRODUCT, B, A; atol=atol)
end

function closure_residuals(ip::AbstractInnerProduct, B::MatrixBasis)
    n = length(B)
    R = Matrix{Float64}(undef, n, n)

    for i in 1:n
        for j in 1:n
            R[i, j] = residual(ip, B, bracket(B, i, j))
        end
    end

    return R
end

function closure_residuals(B::MatrixBasis)
    return closure_residuals(DEFAULT_INNER_PRODUCT, B)
end

function is_closed(ip::AbstractInnerProduct, B::MatrixBasis; atol::Real=1e-10)
    return maximum(closure_residuals(ip, B)) < atol
end

function is_closed(B::MatrixBasis; atol::Real=1e-10)
    return is_closed(DEFAULT_INNER_PRODUCT, B; atol=atol)
end

function assert_closed(ip::AbstractInnerProduct, B::MatrixBasis; atol::Real=1e-10)
    is_closed(ip, B; atol=atol) || throw(ArgumentError("The basis is not closed under the bracket operation."))
    return true
end

function assert_closed(B::MatrixBasis; atol::Real=1e-10)
    return assert_closed(DEFAULT_INNER_PRODUCT, B; atol=atol)
end
