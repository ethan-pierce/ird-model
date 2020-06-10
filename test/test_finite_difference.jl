using Test

include("../src/ModelGrid.jl")
using .ModelGrid: Grid, set!
include("../src/FiniteDifference.jl")
using .FiniteDifference: grad_x, grad_y, lap_x, lap_y, dxdy

@testset "Difference Operators" begin
    g = Grid((10, 10), (1, 1))
    g.data[:, :] .= 1
    set!(g, (x, y) -> sin(x) + cos(y))
    xdiff = grad_x(g)
    println(xdiff.data)
end

