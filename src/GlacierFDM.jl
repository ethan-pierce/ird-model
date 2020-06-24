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
    shape::Tuple{Int, Int}
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

    lastx = grid.shape[1]

    for xj = 1:grid.shape[2]
        output[1, xj] = (1 / (2 * grid.dx)) * (grid[1, xj] - grid[2, xj])
        output[lastx, xj] = (1 / (2 * grid.dx)) * (grid[lastx - 1, xj] - grid[lastx, xj])
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

    lasty = grid.shape[2]

    for xi = 1:grid.shape[1]
        output[xi, 1] = (1 / (2 * grid.dy)) * (grid[xi, 1] - grid[xi, 2])
        output[xi, lasty] = (1 / (2 * grid.dy)) * (grid[xi, lasty - 1] - grid[xi, lasty])
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

function guess_velocity(initial_u::Grid, initial_v::Grid, eta::Grid, H::Grid)::Tuple{Grid, Grid}
    # notation: Lu = f
    # currently only implemented for h = 1
    nx = H.shape[1]
    ny = H.shape[2]
    nn = 2 * (nx * ny)

    U = deepcopy(initial_u)
    V = deepcopy(initial_v)

    Sx = grad_x(H)
    Sy = grad_y(H)
    fx = (rho * g) .* Sx.data
    fy = (rho * g) .* Sy.data
    f = vcat(reshape(fx, :, 1), reshape(fy, :, 1))

    etax = grad_x(eta)
    etay = grad_y(eta)

    L = Array{Float64}(undef, nn, nn)
    
    # Sx component
    # for x = 2:nx - 1
        # for y = 2:ny - 1
            
            
            

    # result = L \ f
    # U = result[1:nn]
    # V = result[nn+1:end]

    return U, V
end

end