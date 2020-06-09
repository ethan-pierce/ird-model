module ModelGrid

mutable struct Grid
    shape::Tuple{UInt32, UInt32, UInt32}
    h::Tuple{UInt32, UInt32, UInt32}
    dx::UInt32
    dy::UInt32
    dz::UInt32
    data::AbstractArray
    Grid(shape, h) = 
        new(shape, h, h[1], h[2], h[3],
            zeros(Float64, shape[1], shape[2], shape[3]))
end


end