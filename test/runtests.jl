import Base.Test: @testset, @test, @test_throws
import StaticArrays: SVector

using MultipleScattering

include("basic_type_tests.jl")
include("shapetests.jl")

include("fourier_tests.jl")

include("acoustic_physics_tests.jl")
include("boundary_conditions_tests.jl")

include("end_to_end_tests.jl")
