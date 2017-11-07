using MultipleScattering
using Base.Test


@testset "Summary" begin

    @testset "Shape" begin

        @testset "Rectangle" begin
            simple_rectangle = Rectangle([0.0,0.0],[2.0,3.0])
            @test volume(simple_rectangle) ≈ 6.0
            @test name(simple_rectangle) == "Rectangle"
        end

        @testset "Circle" begin
            circle = Circle(2.0,[6.7,8.9])
            circle_bounding_box = bounding_box(circle)
            @test volume(circle)/volume(circle_bounding_box) ≈ 0.7853981633974483
            @test name(circle) == "Circle"
        end

        @testset "Time of flight" begin
            time_of_flight = TimeOfFlight([-10.0,0.0],40.0)
            time_of_flight_bounding_box = bounding_box(time_of_flight)
            ratio = volume(time_of_flight)/volume(time_of_flight_bounding_box)
            # Geometric arguments dictate that the ratio must be between 0.5 and 1.0
            @test ratio > 0.5
            @test ratio < 1.0
            @test name(time_of_flight) == "Time of flight"
        end
    end

    @testset "Particle" begin
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

        p1 = Particle([0.0,0.0])
        p1_ptr = p1
        p1b = Particle([-0.0,0.0])

        # isequal is strict about +/-0.0 but IEEE compliant == is not strict
        @test (p1 == p1b)
        @test !isequal(p1,p1b)

        # Check that Julia fallback === and !== work
        @test p1_ptr === p1
        @test p1b !== p1

        # Generate 5 random number
        radii = 1 .+ randn(5).^2
        particles = [ Particle([0.0,0.0],radii[i]) for i=1:5 ]
        @test std_radius(particles) === std(radii)
        @test mean_radius(particles) === mean(radii)
    end

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

    @testset "TimeModel" begin
        # Time response from a single particle
        include("../example/time_model.jl")
        freq_model, time_model = run_time_response_single_particle()
        # Need to test that the spike appears at the right place
        @test true
        
        include("../example/lens.jl")
        freq_model, time_model = run_lens()
        # Need to test that the spike appears at the right place
        @test true
    end

end
