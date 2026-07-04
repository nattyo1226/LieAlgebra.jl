function test_structure_constants()
    B = MatrixBasis([px, py, pz])
    f = structure_constants(B)

    @test f[1, 2, 3] ≈ 2
    @test f[2, 3, 1] ≈ 2
    @test f[3, 1, 2] ≈ 2
    @test is_closed(B)
end

@testset "Structure Constants" begin
    test_structure_constants()
end
