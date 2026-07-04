function test_matrix_basis()
    @test rank(B) == 3
    @test is_independent(B)

    @test in_span(B, commutator(px, py))
    @test bracket_coefficients(B, 1, 2) ≈ [0, 0, 2]
end

@testset "MatrixBasis" begin
    test_matrix_basis()
end
