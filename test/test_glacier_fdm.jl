using Test

include("../src/GlacierFDM.jl")
using .GlacierFDM: Grid, set!, grad_x

@testset "Constructors" begin
    g = Grid((7, 5), (1, 1))
    set!(g, (x, y) -> 1)
    @test size(g.data) == (5, 7)
    @test g.dx == 1
    
    g2 = copy(g)
    @test size(g.data) == (5, 7)
    @test g.dx == 1

    set!(g, (x, y) -> x^2 + y)
    @test g.data[1, 3] == 10
    @test g[3, 1] == 10

    xdiff = grad_x(g)
    @test xdiff[2, 1] == -4.0

end