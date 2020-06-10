using Test

include("../src/ModelGrid.jl")
using .ModelGrid: Grid, set!

@testset "Constructors" begin
    g = Grid((100, 20), (1, 0.5))
    @test length(g.data) == 100 * 20

    f = copy(g)
    @test length(f.data) == 100 * 20

    g.data[:] .= 1
    set!(g, (x, y) -> sin(x) + cos(y))
    @test isapprox(g.data[1, 1], 1.382, atol = 1e-3)
end