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
        # Plot a single particle
        particle = Particle([0.0,0.0])
        plot(particle)

        # This should work, but currently fails
        # # Plot a vector of particles
        # particles = random_particles(10,1.0)
        # plot(particles)
        @test true
    end
end
