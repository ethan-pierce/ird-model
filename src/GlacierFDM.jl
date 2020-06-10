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

Base.copy(grid::Grid) = Grid(grid.shape, grid.h, grid.data)

function set!(grid::Grid, fn::Function)
    for xi = 1:grid.shape[1]
        for xj = 1:grid.shape[2]
            grid[xi, xj] = fn(xi, xj)
        end
    end
    return grid
end

function grad_x(grid::Grid)::Grid
    output = copy(grid)
    for xi = 2:output.shape[1] - 2
        for xj = 1:output.shape[2]
            output[xi, xj] = (1 / (2 * output.dx)) * (output[xi - 1, xj] - output[xi + 1, xj])
        end
    end
    return output
end

function grad_y()
end

function lap_x()
end

function lap_y()
end

function dxdy()
end



end