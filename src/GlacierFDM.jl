module GlacierFDM

mutable struct Grid
    shape::Tuple{UInt32, UInt32}
    h::Tuple{Float32, Float32}
    dx::Float32 
    dy::Float32
    data::Array{Float64, 2}
    Grid(shape, h) = new(shape, h, h[1], h[2], Array{Float64, 2}(undef, shape[2], shape[1]))
    Grid(shape, h, data) = new(shape, h, h[1], h[2], data)
end

function Base.getindex(grid::Grid, xi::Int, xj::Int)
    return grid.data[xj, xi]
end

function Base.setindex!(grid::Grid, val::Real, xi::Int, xj::Int)
    grid.data[xj, xi] = val
    return grid
end

function set!(grid::Grid, fn::Function)
    for xi = 1:grid.shape[1]
        for xj = 1:grid.shape[2]
            grid[xi, xj] = fn(xi, xj)
        end
    end
    return grid
end

function grad_x(grid::Grid)::Grid
    output = deepcopy(grid)
    for xi = 2:output.shape[1] - 1
        for xj = 1:output.shape[2]
            output[xi, xj] = (1 / (2 * output.dx)) * (output[xi - 1, xj] - output[xi + 1, xj])
        end
    end
    return output
end

function grad_y(grid::Grid)::Grid
    output = deepcopy(grid)
    for xi = 1:output.shape[1]
        for xj = 2:output.shape[2] - 1
            output[xi, xj] = (1 / (2 * output.dy)) * (output[xi, xj - 1] - output[xi, xj + 1])
        end
    end
    return output
end

function calculate_strain(u::Grid, v::Grid)::Grid
    epsilon = deepcopy(u)
    dudy = grad_y(u)
    dvdx = grad_x(v)
    epsilon.data = (1 / 2) .* (dudy.data + dvdx.data)
    return epsilon
end

function B(temperature::Grid, H::Grid)::Grid
    const rho = 917 # kg m^-3
    const g = 9.81 # m s^-2
    const phi = 9.8e-8 # K Pa^-1
    const A_0 = 2.4e-24 # s Pa^-n 
    const n = 3 
    const Q = 115e3 # J mol^-1
    const R = 8.3145 # J mol^-1 K^-1
    T_adj = temperature.data .+ (rho * g * phi) .* H.data
    A_array = A_0 .* exp(-Q ./ (R .* T_adj))
    B_array = A_array.^(-1 / n)
    B = deepcopy(temperature)
    B.data = B_array
    return B
end

function calculate_eta(strain::Grid, temperature::Grid, H::Grid)::Grid
    const n = 3
    B = B(temperature, H)
    strain_norm = (1 / 2) .* strain.data .* strain.data
    eta = deepcopy(temperature)
    eta.data = (1 / 2) .* B.data .* (strain_norm .^ ((1 - n) / n))
end

end