using MultipleScattering
using Base.Test


@testset "Summary" begin
    # @test single_scatter_test()

    @testset "Shape tests" begin
        # Test shapes and their bounding boxes
        circle = Circle(2.0,[6.7,8.9])
        circle_bounding_box = bounding_box(circle)

        @test volume(circle)/volume(circle_bounding_box) ≈ 0.7853981633974483

        time_of_flight = TimeOfFlight([-10.0,0.0],40.0)
        time_of_flight_bounding_box = bounding_box(time_of_flight)
        ratio = volume(time_of_flight)/volume(time_of_flight_bounding_box)

        # Geometric arguments dictate that the ratio must be between 0.5 and 1.0
        @test ratio > 0.5
        @test ratio < 1.0

        simple_rectangle = Rectangle([0.0,0.0],[2.0,3.0])
        @test volume(simple_rectangle) ≈ 6.0
    end

    @testset "Particle test" begin
        # Make two random seeds, extremely low probability they will be the same
        seed1 = Base.Random.make_seed()
        seed2 = Base.Random.make_seed()

        volfrac = 0.2
        radius = 0.5
        shape = Circle(10.0,[0.0,0.0])

        particles1 = random_particles(volfrac, radius, shape; seed = seed1)
        particles1a = random_particles(volfrac, radius, shape; seed = seed1)
        particles2 = random_particles(volfrac, radius, shape; seed = seed2)

        # Particles should be determined solely by the seed
        @test particles1 == particles1a
        @test particles1 != particles2
    end

    @testset "Scattering test" begin
        include("single_scatter.jl")
        # Test against analytic solution
        @test single_scatter_test()
    end

    @testset "Plot test" begin
        using MultipleScattering.Plot
        # Just run it to see if we have any errors (yes thats a very low bar)

        volfrac = 0.01
        radius = 2.0
        k_arr = collect(linspace(0.2,1.0,5))
        model = FrequencyModel(volfrac,radius,k_arr)

        plot_model(model)
        plot_field(model,0.2)

        @test true
    end

end
