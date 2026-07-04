function test_decompose_su2()
    B = MatrixBasis([px, py, pz])

    rng = MersenneTwister(816)
    D = decompose(rng, B)

    @test length(D.components) == 1
    @test sort(D.components[1]) == [1, 2, 3]

    B2 = transformed_basis(B, D)

    @test length(B2) == 3
    @test is_closed(B2)
end

function test_decompose_su2_plus_su2()
    Sx = kron(px, I2)
    Sy = kron(py, I2)
    Sz = kron(pz, I2)

    Tx = kron(I2, px)
    Ty = kron(I2, py)
    Tz = kron(I2, pz)

    B = MatrixBasis([Sx, Sy, Sz, Tx, Ty, Tz])

    rng = MersenneTwister(816)
    D = decompose(rng, B)

    @test length(D.components) == 2
    @test sort(length.(D.components)) == [3, 3]

    B2 = transformed_basis(B, D)
    comps = component_bases(B2, D)

    @test length(B2) == 6
    @test length(comps) == 2
    @test all(length.(comps) .== 3)

    @test all(C -> is_closed(C), comps)

    C1, C2 = comps

    for A in C1
        for B in C2
            @test norm(lie_bracket(A, B)) < 1e-10
        end
    end
end

@testset "Decomposition" begin
    test_decompose_su2()
    test_decompose_su2_plus_su2()
end
