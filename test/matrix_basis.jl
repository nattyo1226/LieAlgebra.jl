function test_matrix_basis()
    B = MatrixBasis([px, py, pz])

    @test rank(B) == 3
    @test is_independent(B)

    @test in_span(B, commutator(px, py))
    @test bracket_coefficients(B, 1, 2) ≈ [0, 0, 2im]

    f = structure_constants(B)
    @test f[1, 2, 3] ≈ 2im
    @test f[2, 3, 1] ≈ 2im
    @test f[3, 1, 2] ≈ 2im
    @test is_closed(B)
end

@testset "MatrixBasis" begin
    test_matrix_basis()
end
