module FiniteDifference

include("./ModelGrid.jl")
using .ModelGrid: Grid 

function dx_1(grid::Grid)::AbstractArray
    result = zeros(grid)
    for y = 2:grid.shape[2] - 1
        for x = 2:grid.shape[1] - 1
            for z = 2:grid.shape[3] - 1
                a = 0
            end
        end
    end
end








end