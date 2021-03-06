using Test

include("../src/GlacierFDM.jl")
using .GlacierFDM: Grid, set!, fill_zeros!, grad_x, grad_y,
    calculate_strain, B, calculate_eta, guess_velocity

@testset "Constructors" begin
    g = Grid((7, 5), (1, 1))
    set!(g, (x, y) -> 1)
    @test size(g.data) == (5, 7)
    @test g.dx == 1
    
    g2 = deepcopy(g)
    @test size(g2.data) == (5, 7)
    @test g2.dx == 1

    g[1, 1] = 100
    @test g2[1, 1] != g[1, 1] 

    set!(g, (x, y) -> x^2 + y)
    @test g.data[1, 3] == 10
    @test g[3, 1] == 10
end

@testset "Gradients" begin
    g = Grid((7, 5), (1, 1))
    set!(g, (x, y) -> x^2 + y)

    xdiff = grad_x(g)
    @test xdiff[2, 1] == -4.0
    @test xdiff[1, 4] == -1.5

    ydiff = grad_y(g)
    @test ydiff[1, 2] == -1
    @test ydiff[3, 1] == -0.5
end

@testset "Effective viscosity" begin
    u = Grid((4, 4), (1, 1))
    set!(u, (x, y) -> 3*(x + y))

    v = Grid((4,4), (1, 1))
    set!(v, (x, y) -> 1 * (x + y))
    
    epsilon = calculate_strain(u, v)
    @test epsilon[2, 2] == -2.0
    @test epsilon[2, 2] == epsilon[2, 3]
    @test epsilon[2, 2] == epsilon[3, 2]

    T = Grid((4, 4), (1, 1))
    set!(T, (x, y) -> 265 + y + x)

    H = Grid((4, 4), (1, 1))
    set!(H, (x, y) -> 1000 + 2 * x)

    Btest = B(T, H)
    eta = calculate_eta(epsilon, T, H)
end

@testset "Velocity snapshot" begin
    u0 = Grid((4, 4), (1, 1))
    v0 = Grid((4, 4), (1, 1))
    T = Grid((4, 4), (1, 1))
    H = Grid((4, 4), (1, 1))

    set!(u0, (x, y) -> x)
    set!(v0, (x, y) -> 1e-2 * x)
    set!(T, (x, y) -> 263 + y + x)
    set!(H, (x, y) -> 1e3 - 2 * x)
    
    epsilon = calculate_strain(u0, v0)
    eta = calculate_eta(epsilon, T, H)
    u, v = guess_velocity(u0, v0, eta, H)
    println(u)
end