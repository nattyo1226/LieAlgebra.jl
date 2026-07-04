function test_structure_constants()
    f = structure_constants(B)

    @test f[1, 2, 3] ≈ 2im
    @test f[2, 3, 1] ≈ 2im
    @test f[3, 1, 2] ≈ 2im
    @test is_closed(B)
end

@testset "Structure Constants" begin
    test_structure_constants()
end
