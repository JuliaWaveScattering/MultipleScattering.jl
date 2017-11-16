using MultipleScattering
using Base.Test
using Plots

@testset "Summary" begin

    include("shapetests.jl")
    include("particle_tests.jl")
    include("time_model_tests.jl")
    include("moments_tests.jl")

    @testset "Type stability" begin
        # Define everything as a Float32
        volfrac = 0.01f0
        radius = 1.0f0
        k_arr = collect(linspace(0.01f0,1.0f0,100))
        model = FrequencyModel(volfrac,radius,k_arr)
        @test typeof(model.response[1]) == Complex{Float32}
    end

    @testset "Single scatterer" begin
        include("single_scatter.jl")
        # Test against analytic solution
        @test single_scatter_test()
    end

    @testset "Boundary conditions" begin
        include("boundary_conditions.jl")
        # Test boundary conditions for 4 particles with random properties and positions.
        @test boundary_conditions_test()
    end

    @testset "Plot FrequencyModel" begin
        # Just run it to see if we have any errors (yes thats a very low bar)
        volfrac = 0.01
        radius = 2.0
        k_arr = collect(linspace(0.2,1.0,5))
        model = FrequencyModel(volfrac,radius,k_arr)

        plot(model)
        plot(model,0.2)

        @test true
    end

end
