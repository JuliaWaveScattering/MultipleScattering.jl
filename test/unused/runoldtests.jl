using MultipleScattering
using Base.Test
using Plots

@testset "Summary" begin

    include("shapetests.jl")
    include("particle_tests.jl")
    include("time_simulation_tests.jl")
    include("moments_tests.jl")

    # test our discretised Fourier transforms
    @testset "Fourier tranforms" begin
      include("fourier_test.jl")
      @test fourier_test()
    end

    @testset "Type stability" begin
        # Define everything as a Float32
        volfrac = 0.01f0
        radius = 1.0f0
        k_arr = collect(LinRange(0.01f0,1.0f0,100))
        simulation = FrequencySimulation(volfrac,radius,k_arr)
        @test typeof(simulation.response[1]) == Complex{Float32}
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

    @testset "Plot FrequencySimulation" begin
        # Just run it to see if we have any errors (yes thats a very low bar)
        volfrac = 0.01
        radius = 2.0
        k_arr = collect(LinRange(0.2,1.0,5))
        simulation = FrequencySimulation(volfrac,radius,k_arr)

        plot(simulation)
        plot(simulation,0.2)

        @test true
    end

end
