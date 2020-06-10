using SafeTestsets

# @safetestset "Model Grid" begin include("test_model_grid.jl") end
# @safetestset "Finite Difference" begin include("test_finite_difference.jl") end

@safetestset "Glacier FDM" begin include("test_glacier_fdm.jl") end