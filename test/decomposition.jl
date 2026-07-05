function test_decompose_su2()
    B = MatrixBasis([px, py, pz])

    rng = MersenneTwister(816)
    D = decompose(rng, B)

    @test length(D.ideals) == 1
    @test length(D[1]) == 3
    @test is_closed(D[1])
    @test all(ishermitian, D[1])
end

function test_decompose_su2_plus_su2()
    B = MatrixBasis([Sx, Sy, Sz, Tx, Ty, Tz])

    rng = MersenneTwister(816)
    D = decompose(rng, B)

    @test length(D) == 2
    @test sort(length.(D)) == [3, 3]
    @test all(is_closed, D)
    @test all(all(ishermitian, ideal) for ideal in D)

    for B1 in D[1]
        for B2 in D[2]
            @test norm(lie_bracket(B1, B2)) < 1e-10
        end
    end
end

@testset "Decomposition" begin
    test_decompose_su2()
    test_decompose_su2_plus_su2()
end
