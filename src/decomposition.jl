Base.length(D::Decomposition) = length(D.ideals)
Base.getindex(D::Decomposition, i::Integer) = D.ideals[i]
Base.iterate(D::Decomposition, state...) = iterate(D.ideals, state...)
Base.keys(D::Decomposition) = Base.OneTo(length(D))
Base.eachindex(D::Decomposition) = Base.OneTo(length(D))

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

function transform_basis(B::MatrixBasis, U::AbstractMatrix)
    n = length(B)
    size(U) == (n, n) || throw(DimensionMismatch("Transformation matrix must be $(n) × $(n)."))

    elements = [
        reconstruct(B, U[:, i])
        for i in 1:n
    ]

    return MatrixBasis(elements)
end

function component_basis(B::MatrixBasis, inds::AbstractVector{<:Integer})
    return MatrixBasis([B[i] for i in inds])
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
    B2 = transform_basis(B, F.vectors)

    ideals = [component_basis(B2, group) for group in groups]

    return Decomposition(ideals)
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
