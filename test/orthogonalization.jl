function test_orthogonalization_1()
    B = MatrixBasis([px, py, pz])
    B2 = orthogonalize(B)

    @test length(B2) == 3
    @test is_orthonormal(B2)
    @test is_closed(B2)
    @test all(ishermitian, B2)
end

function test_orthogonalization_2()
    B = MatrixBasis([Sx, Sy, Sz, Tx, Ty, Tz])

    rng = MersenneTwister(816)
    D = decompose(rng, B)
    D2 = orthogonalize(D)

    @test length(D2) == 2
    @test all(is_orthonormal, D2)
    @test all(is_closed, D2)
    @test all(all(ishermitian, ideal) for ideal in D2)
end

@testset "Orthogonalization" begin
    test_orthogonalization_1()
    test_orthogonalization_2()
end
