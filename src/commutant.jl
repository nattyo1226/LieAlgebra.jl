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

function commutant_complex(As::AbstractVector{<:AbstractMatrix}; atol::Real=1e-10)
    M = commutant_matrix(As)
    ns = nullspace(M; atol=atol)

    n = size(first(As), 1)

    return [reshape(ns[:, i], (n, n)) for i in axes(ns, 2)]
end

function real_matrix(A::AbstractMatrix; atol::Real=1e-10)
    max_im = maximum(abs.(imag.(A)))
    max_im < atol || throw(ArgumentError("Matrix is not real within tolerance $(atol)."))
    return real.(A)
end

function hermitian_part(A::AbstractMatrix)
    return (A + A') / 2
end

function antihermitian_part_to_hermitian(A::AbstractMatrix)
    return (A - A') / (2im)
end

function vec_real(A::AbstractMatrix)
    a = vec(A)
    return vcat(real(a), imag(a))
end

function independent_matrices_real(As::AbstractVector{<:AbstractMatrix}; atol::Real=1e-10)
    isempty(As) && return Vector{Matrix{Float64}}()

    V = hcat(vec_real.(As)...)
    F = qr(V, ColumnNorm())

    r = sum(abs.(diag(F.R)) .> atol)
    inds = F.p[1:r]

    return [Matrix{Float64}(As[i]) for i in inds]
end

function commutant(As::AbstractVector{<:AbstractMatrix}; atol::Real=1e-10)
    Cs = commutant_complex(As; atol=atol)
    Hs = Vector{Matrix{Float64}}()

    for C in Cs
        h1 = hermitian_part(C)
        h2 = antihermitian_part_to_hermitian(C)

        push!(Hs, real_matrix(h1; atol=atol))
        push!(Hs, real_matrix(h2; atol=atol))
    end

    return independent_matrices_real(Hs; atol=atol)
end

function commutant(ip::AbstractInnerProduct, B::MatrixBasis; atol::Real=1e-10)
    As = adjoint_matrices(ip, B)
    return commutant(As; atol=atol)
end

function commutant(B::MatrixBasis; atol::Real=1e-10)
    return commutant(DEFAULT_INNER_PRODUCT, B; atol=atol)
end
