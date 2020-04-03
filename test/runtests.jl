import Test: @testset, @test, @test_throws
import StaticArrays: SVector

using MultipleScattering
using LinearAlgebra
using Plots

include("special_functions.jl")

include("shapetests.jl")
include("particle_tests.jl")

include("print_tests.jl")

include("fourier_tests.jl")

include("acoustic_physics_tests.jl")
include("source_test.jl")
include("boundary_conditions_tests.jl")

include("end_to_end_tests.jl")
include("time_tests.jl")

include("moments_tests.jl")

include("run_examples.jl")
