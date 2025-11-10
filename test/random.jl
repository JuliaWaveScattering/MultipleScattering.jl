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


@testset "Box around particles" begin
    seed1 = 1

    volfrac = 0.1
    radius = 0.5

    particle_shape = Sphere(radius)
    medium = Acoustic(3; ρ = 0.2, c = 0.2)

    region_shape = Box([0.0,0.0,0.0],[10.0,12.0,7.0])

    particles = random_particles(medium, particle_shape, region_shape, volfrac; seed=seed1)
    box = bounding_box(particles)

    # check that all particles are within the bounding box 
    @test all(p ⊆ box for p in particles)

    # check that moving all particles makes at least one particle leave the box
    particles2 = [
        Particle(medium, 
            congruent(particle_shape, origin(p) + 100 .* [eps(Float64),eps(Float64), eps(Float64)])
        )
    for p in particles]    

    @test any(!(p ⊆ box) for p in particles2)

    # create a box around particle centres
    points = origin.(particles)
    box = Box(points)

    @test all(p ∈ box for p in points)
    @test !all(p ⊆ box for p in particles)

    points = [
        p + 100 .* [eps(Float64), eps(Float64), eps(Float64)] 
    for p in points]
    @test !all(p ∈ box for p in points)
end