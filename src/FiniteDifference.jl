module FiniteDifference

include("./ModelGrid.jl")
using .ModelGrid: Grid

function grad_x(grid::Grid)::Grid
    output = copy(grid)
    for y = 2:output.shape[2] - 2
        for x = 2:output.shape[1] - 2
            output[x, y] = (1 / (2 * output.h[1])) * (output[x - 1, y] - output[x + 1, y])
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