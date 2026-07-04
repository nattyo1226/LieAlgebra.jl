function generic_commutant_element(rng::AbstractRNG, Cs::AbstractVector{<:AbstractMatrix})
    _check_empty(Cs)
    _check_same_square(Cs)

    X = zero(first(Cs))

    for C in Cs
        X += randn(rng) * C
    end

    return Hermitian(X)
end

function generic_commutant_element(Cs::AbstractVector{<:AbstractMatrix})
    return generic_commutant_element(Random.default_rng(), Cs)
end

function eigenspace_groups(vals::AbstractVector; atol::Real=1e-10)
    groups = Vector{Vector{Int}}()

    for i in eachindex(vals)
        found = false

        for group in groups
            if abs(vals[i] - vals[first(group)]) < atol
                push!(group, i)
                found = true
                break
            end
        end

        found || push!(groups, [i])
    end

    return groups
end

function decompose(
    rng::AbstractRNG,
    ip::AbstractInnerProduct,
    B::MatrixBasis;
    atol::Real=1e-10,
)
    Cs = commutant(ip, B; atol=atol)
    X = generic_commutant_element(rng, Cs)

    F = eigen(X)
    groups = eigenspace_groups(F.values; atol=atol)

    return Decomposition(groups, Matrix{ComplexF64}(F.vectors))
end

function decompose(
    ip::AbstractInnerProduct,
    B::MatrixBasis;
    atol::Real=1e-10,
)
    return decompose(Random.default_rng(), ip, B; atol=atol)
end

function decompose(
    rng::AbstractRNG,
    B::MatrixBasis;
    atol::Real=1e-10,
)
    return decompose(rng, DEFAULT_INNER_PRODUCT, B; atol=atol)
end

function decompose(
    B::MatrixBasis;
    atol::Real=1e-10,
)
    return decompose(Random.default_rng(), DEFAULT_INNER_PRODUCT, B; atol=atol)
end

function transform_basis(B::MatrixBasis, U::AbstractMatrix)
    n = length(B)
    size(U) == (n, n) || throw(DimensionMismatch("Transformation matrix must be $(n) × $(n)."))

    new_elements = [
        reconstruct(B, U[:, i])
        for i in 1:n
    ]

    return MatrixBasis(new_elements)
end

function transformed_basis(B::MatrixBasis, D::Decomposition)
    return transform_basis(B, D.transform)
end

function component_basis(B::MatrixBasis, component::AbstractVector{<:Integer})
    return MatrixBasis([B[i] for i in component])
end

function component_bases(B::MatrixBasis, D::Decomposition)
    B2 = transformed_basis(B, D)
    return [component_basis(B2, component) for component in D.components]
end
