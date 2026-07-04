function test_commutant()
    B = MatrixBasis([px, py, pz])
    C = commutant(B)

    @test length(C) == 1
    @test C[1] ≈ Matrix(I(3) / sqrt(3))
end

@testset "Commutant" begin
    test_commutant()
end
