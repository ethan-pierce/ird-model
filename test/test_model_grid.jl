using Test

include("../src/ModelGrid.jl")
using .ModelGrid: Grid

@testset "Constructor" begin
    g = Grid((100, 20, 10), (10, 2, 1))
    @test length(g.data) == 100 * 20 * 10

    g.data[1, 1, 1] = 1
    @test g.data[1, 1, 1] == 1
end