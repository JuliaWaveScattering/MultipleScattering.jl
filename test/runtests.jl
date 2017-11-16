using MultipleScattering
using Base.Test
using Plots

@testset "Summary" begin

    @testset "Shape" begin

        @testset "Rectangle" begin
            rectangle = Rectangle([0.0,0.0],[2.0,3.0])
            @test volume(rectangle) ≈ 6.0
            @test name(rectangle) == "Rectangle"

            @testset "Boundary functions" begin
                x, y = boundary_functions(rectangle)
                @test x.(0:0.1:1) ≈ [0.0, 0.8, 1.6, 2.0, 2.0, 2.0, 1.2, 0.4, 0.0, 0.0, 0.0]
                @test y.(0:0.1:1) ≈ [0.0, 0.0, 0.0, 0.6, 1.8, 3.0, 3.0, 3.0, 2.4, 1.2, 0.0]
                @test_throws(Exception,x(-eps(Float64)))
                @test_throws(Exception,x(1.0 + eps(Float64)))
                @test_throws(Exception,y(-eps(Float64)))
                @test_throws(Exception,y(1.0 + eps(Float64)))
            end
        end

        @testset "Circle" begin
            circle = Circle(2.0,[6.7,8.9])
            circle_bounding_box = bounding_box(circle)
            @test volume(circle)/volume(circle_bounding_box) ≈ 0.7853981633974483
            @test name(circle) == "Circle"

            @testset "Boundary functions" begin
                x, y = boundary_functions(Circle(3.0,[-1.0,2.0]))
                @test x.(0:0.1:1) ≈ [0.0,-0.19098300562505255,-0.6909830056250525,-1.3090169943749475,-1.8090169943749475,-2.0,-1.8090169943749475,-1.3090169943749475,-0.6909830056250528,-0.19098300562505266,0.0]
                @test y.(0:0.1:1) ≈ [2.0,2.5877852522924734,2.9510565162951536,2.9510565162951536,2.5877852522924734,2.0,1.412214747707527,1.0489434837048464,1.0489434837048464,1.4122147477075266,1.9999999999999998]
                @test_throws(Exception,x(-eps(Float64)))
                @test_throws(Exception,x(1.0 + eps(Float64)))
                @test_throws(Exception,y(-eps(Float64)))
                @test_throws(Exception,y(1.0 + eps(Float64)))
            end
        end

        @testset "Time of flight" begin
            time_of_flight = TimeOfFlight([-10.0,0.0],40.0)
            time_of_flight_bounding_box = bounding_box(time_of_flight)
            ratio = volume(time_of_flight)/volume(time_of_flight_bounding_box)
            # Geometric arguments dictate that the ratio must be between 0.5 and 1.0
            @test ratio > 0.5
            @test ratio < 1.0
            @test name(time_of_flight) == "Time of flight"

            @testset "Boundary functions" begin
                x, y = boundary_functions(TimeOfFlight([-1.0,0.0],3.0))
                @test x.(0:0.1:1) ≈ [0.0,1.2182846930812632,1.9095426112571139,1.9095426112571139,1.2182846930812632,0.0,0.0,0.0,0.0,0.0,0.0]
                @test y.(0:0.1:1) ≈ [-2.82842712474619,-2.019706171808505,-0.7311373286046446,0.7311373286046446,2.0197061718085054,2.82842712474619,1.6970562748477145,0.5656854249492386,-0.5656854249492386,-1.6970562748477145,-2.8284271247461903]
                @test_throws(Exception,x(-eps(Float64)))
                @test_throws(Exception,x(1.0 + eps(Float64)))
                @test_throws(Exception,y(-eps(Float64)))
                @test_throws(Exception,y(1.0 + eps(Float64)))
            end
        end

        @testset "Plot Shapes" begin
            # Just run it to see if we have any errors (yes thats a very low bar)
            include("../example/plot_shapes.jl")
            @test true
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

        @testset "Plot Particles" begin
            # Just run it to see if we have any errors (yes thats a very low bar)
            include("../example/plot_particles.jl")
            @test true
        end
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
        # Just run it to see if we have any errors (yes thats a very low bar)
        volfrac = 0.01
        radius = 2.0
        k_arr = collect(linspace(0.2,1.0,5))
        model = FrequencyModel(volfrac,radius,k_arr)

        plot(model)
        plot(model,0.2)

        @test true
    end

    @testset "TimeModel" begin
        # Time response from a single particle
        include("../example/time_model.jl")
        freq_model, time_model = run_time_response_single_particle()
        # Spike at start and then two diffracted spikes later, then
        # not much else 

              # incident plane wave
        @test abs(time_model.response[1] - 16) < 1e-1     &&
              # first reflection
              time_model.response[3473] > 0.7             &&
              abs(time_model.time_arr[3473] - 218) < 1e-1 &&
              # largest diffraction
              time_model.response[3601] > 0.9             &&
              abs(time_model.time_arr[3601] - 226) < 1e-1

        # Take samples every 200 and compare to
        # previously run result
        @test time_model.response[1:200:4000,1] ≈ [15.91568766260155,-0.010114700804103129,-0.010849657057909082,-0.0094274294391117,-0.008064690305764446,-0.008065776499280565,-0.008227945663215533,-0.007811684591217631,-0.0066031417092440635,-0.005390619379328917,-0.004513676793355925,-0.003388964946642101,-0.0017444062810255325,-9.962160860650897e-5,0.0011156638442981646,0.002269075560701978,0.003948284351502541,0.007065202861450666,0.9172645743771086,0.014864642396833697]

        include("../example/lens.jl")
        freq_model, time_model = run_lens()
        # Need to test that the spike appears at the right place
        @test true
    end

    @testset "Moments" begin
        begin
            # Test against a problem which can be easily solved
            models = Vector{FrequencyModel{Float64}}(3)
            # Fake responses, with mean 4.0, standard deviation 2.0 and skew 0.0
            responses = [2.0, 4.0, 6.0]
            particles = [Particle([0.0,0.0])]
            for i=1:3
                models[i] = FrequencyModel(particles,[1.0];generate_responses=false)
                models[i].response = reshape([responses[i]+0.0im],1,1)
            end
            moments = Moments(models)
            @test moments.moments[1][1] ≈ 4.0 &&
                  moments.moments[2][1] ≈ 2.0 &&
                  moments.moments[3][1] ≈ 0.0
        end
        begin
            # Test against a previously computed problem with a known seed
            include("../example/moments.jl")
            moments = moments_example()
            @test moments.moments[1][23] ≈ 0.9323251967911877  &&
                  moments.moments[2][78] ≈ 0.2487558908409563  &&
                  moments.moments[3][91] ≈ 0.15597075451712927 &&
                  moments.moments[4][32] ≈ 0.17865717302214346
        end
    end

end
