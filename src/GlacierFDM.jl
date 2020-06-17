module GlacierFDM

const zero = 1e-30
const rho = 917 # kg m^-3
const g = 9.81 # m s^-2
const phi_CISM = 9.8e-8 # K Pa^-1
const phi_CP = 7e-8 # K Pa^-1
const A_0 = 2.4e-24 # s Pa^-n 
const n = 3 
const Q_plus = -1.15e5 # J mol^-1
const Q_minus = -6e4 # J mol^-1
const R = 8.3145 # J mol^-1 K^-1
const transition_temperature = 263 # K

mutable struct Grid
    shape::Tuple{UInt32, UInt32}
    h::Tuple{Float32, Float32}
    dx::Float32 
    dy::Float32
    data::Array{BigFloat, 2}
    Grid(shape, h) = new(shape, h, h[1], h[2], Array{BigFloat, 2}(undef, shape[2], shape[1]))
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

function fill_zeros!(grid::Grid)
    for i = 1:length(grid.data)
        if grid.data[i] == 0
            grid.data[i] = zero
        end
    end
    return grid
end

function grad_x(grid::Grid)::Grid
    output = deepcopy(grid)
    for xi = 2:grid.shape[1] - 1
        for xj = 1:grid.shape[2]
            output[xi, xj] = (1 / (2 * grid.dx)) * (grid[xi - 1, xj] - grid[xi + 1, xj])
        end
    end
    return output
end

function grad_y(grid::Grid)::Grid
    output = deepcopy(grid)
    for xi = 1:grid.shape[1]
        for xj = 2:grid.shape[2] - 1
            output[xi, xj] = (1 / (2 * grid.dy)) * (grid[xi, xj - 1] - grid[xi, xj + 1])
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

function set_creep_activation(T_effective::Array, T_transition::Array)::Array
    if length(T_effective) != length(T_transition) && 
        throw(DimensionMismatch("temperature arrays are not the same length"))
    end
    
    Q_array = deepcopy(T_effective)

    for i = 1:length(Q_array)
        if T_effective[i] > T_transition[i]
            Q_array[i] = Q_plus
        else
            Q_array[i] = Q_minus
        end
    end

    return Q_array
end

function B(temperature::Grid, H::Grid)::Grid 
    T_effective = temperature.data .+ (rho * g * phi_CP) .* H.data
    T_transition = H.data .* (rho * g * phi_CP) .+ transition_temperature

    Q_array = set_creep_activation(T_effective, T_transition)
    A_array = A_0 .* exp.(Q_array ./ (R .* T_effective))

    B_array = A_array .^ (-1 / n)
    B = deepcopy(temperature)
    B.data = B_array
    return B
end

function calculate_eta(strain::Grid, temperature::Grid, H::Grid)::Grid
    B_array = B(temperature, H)
    strain_norm = (1 / 2) .* strain.data .* strain.data

    for i = 1:length(strain_norm)
        if strain_norm[i] == 0
            strain_norm[i] = zero
        end
    end

    eta = deepcopy(H)
    eta.data = (1 / 2) .* B_array.data .* (strain_norm .^ ((1 - n) / n))
    return eta
end

end