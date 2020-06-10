module ModelGrid

mutable struct Grid
    shape::Tuple{UInt32, UInt32}
    h::Tuple{Float32, Float32}
    dx::Float32 
    dy::Float32
    data::Array{Float64, 2}
    Grid(shape, h) = new(shape, h, h[1], h[2], Array{Float64, 2}(undef, shape[1], shape[2]))
end

Base.copy(grid::Grid) = Grid(grid.shape, grid.h)

function set!(grid::Grid, fn::Function)
    for y = 1:grid.shape[2]
        for x = 1:grid.shape[1]
            grid.data[x, y] = fn(x, y)
        end
    end
    return grid
end

end