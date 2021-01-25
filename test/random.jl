@testset "Random generation" begin

    # Make two random seeds, extremely low probability they will be the same
    seed1 = 1
    seed2 = 2

    volfrac = 0.1
    radius = 0.5

    region_shape = Sphere(2, 20.0)
    particle_shape = Circle(radius)
    medium = Acoustic(2; ρ = 0.2, c = 0.2)

    particles1  = random_particles(medium, particle_shape, region_shape, volfrac; seed=seed1)
    particles1a = random_particles(medium, particle_shape, region_shape, volfrac; seed=seed1)
    particles2  = random_particles(medium, particle_shape, region_shape, volfrac; seed=seed2)

    # Particles should be determined solely by the seed
    @test particles1 == particles1a
    @test particles1 != particles2

    @test_throws ErrorException random_particles(medium, particle_shape, region_shape, 0.9)
    @test_throws ErrorException random_particles(medium, particle_shape, region_shape, 0.5; seed=3)

    region_shape = Sphere(3, 20.0)
    particle_shape = Sphere(radius)
    medium = Acoustic(3; ρ = 0.2, c = 0.2)

    # If a seed is not provided, then a different seed will be used everytime
    particles2  = random_particles(medium, particle_shape, region_shape, volfrac)
end
