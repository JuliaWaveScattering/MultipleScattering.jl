@testset "Particle" begin

    a = Acoustic(1.0, 1.0, 2)
    a2 = Acoustic(1.0, 2.0, 2)

    @testset "Types" begin

        circle1 = Circle((0.0,0.0),1.0)
        circle2 = Circle((0.0,5.0),2.0)
        a = Acoustic(1.0,1.0,2)
        homog_particles = [Particle(a,circle1), Particle(a,circle2)]

        # Check types comparisons work as user would expect
        @test typeof(homog_particles) <: AbstractParticles
        @test typeof(homog_particles) <: AbstractParticles{Float64}
        @test typeof(homog_particles) <: AbstractParticles{Float64,2}

        circle = Circle((0.0,0.0),1.0)
        rect = Rectangle((2.0,2.0),3.0,2.0)
        diff_shape_particles = [Particle(a,circle), Particle(a,rect)]

        # Check types comparisons work as user would expect
        @test typeof(diff_shape_particles) <: AbstractParticles
        @test typeof(diff_shape_particles) <: AbstractParticles{Float64}
        @test typeof(diff_shape_particles) <: AbstractParticles{Float64,2}

        a2 = Acoustic(1.0,1.0,2)
        a3 = Acoustic(1.0,1.0,3)
        sphere = Sphere((0.0,0.0,0.0),1.0)

        circular_particle = Particle(a2,circle)
        @test typeof(circular_particle) <: Particle{Float64,2}

        spherical_particle = Particle(a3,sphere)
        @test typeof(spherical_particle) <: AbstractParticle{Float64,3}

        # Dimension mismatch throws error
        @test_throws MethodError Particle(a3,circle)
        @test_throws MethodError Particle(a2,sphere)

        # This is a valid vector of valid particles, but not of type Particles
        # because the dimensions don't match
        invalid_particles = [circular_particle, spherical_particle]
        @test_throws TypeError invalid_particles::(Vector{Pt} where {Dim, Pt <: AbstractParticle{Float64,Dim}})

        # does not throw an error
        diff_shape_particles::(Vector{Pt} where {Dim, Pt <: AbstractParticle{Float64,Dim}})
        @test true

    end


    @testset "Comparisons" begin

        circle = Circle((1.0,3.0),2.0)
        circle_identical = Circle((1.0,3.0),2.0)
        circle_congruent = Circle((4.0,7.0),2.0)
        rect = Rectangle((1.0,2.0),2.0,4.0)

        # Construct three particles, with two the same
        p = Particle(a2,circle)
        p_reference = p
        p_identical = Particle(a2,circle_identical)
        p_different = Particle(a2,rect)
        p_congruent = Particle(a2,circle_congruent)

        # Test comparison operators
        @test p == p_identical
        @test p != p_different
        @test !(p == p_different)
        @test iscongruent(p, p_congruent)
        @test !iscongruent(p, p_different)

        # Check that Julia fallback === works
        @test p === p_reference
        # Particles are immutable, so are compared at the bit levelS
        @test p === p_identical

        # Construct two almost identical particles
        p1 = Particle(a, Circle((0.0,0.0),1.0))
        p1b = Particle(a, Circle((-0.0,0.0),1.0))

        # isequal is strict about +/-0.0 but IEEE compliant == is not strict
        @test p1 == p1b
        @test !isequal(p1,p1b)

    end


    @testset "Plot" begin
        # Just run it to see if we have any errors (yes thats a very low bar)

        # Single particle
        particle = Particle(a, Circle((1.0,0.0), 2.0))
        plot(particle)

        # Vector of particles of same shape
        particles = [
            Particle(a, Circle((1.0,0.0), 2.0)),
            Particle(a, Circle((3.0,-1.0), 1.0))
        ]
        plot(particles)

        # Vector of particles of different shape
        particles = [
            CapsuleParticle(
                Particle(Acoustic(2.0,2.0,2),Circle((5.0,3.0),1.0)),
                Particle(Acoustic(1.0,1.,2),Circle((5.0,3.0),2.0))
            ),
            Particle(a, Rectangle((1.0,0.0), 2.0, 3.0)),
            Particle(a, Circle((3.0,-1.0), 1.0))
        ]
        plot(particles)

        @test true
    end


    @testset "Random generation" begin

        # Make two random seeds, extremely low probability they will be the same
        seed1 = 1
        seed2 = 2

        region_shape = Circle(20.0)

        volfrac = 0.1
        radius = 0.5
        particle_shape = Circle(radius)
        medium = Acoustic(1.0,1.0,2)
        particles1  = random_particles(medium, particle_shape, region_shape, volfrac; seed=seed1)
        particles1a = random_particles(medium, particle_shape, region_shape, volfrac; seed=seed1)
        particles2  = random_particles(medium, particle_shape, region_shape, volfrac; seed=seed2)

        # Particles should be determined solely by the seed
        @test particles1 == particles1a
        @test particles1 != particles2

        @test_throws ErrorException random_particles(medium, particle_shape, region_shape, 0.9)
        @test_throws ErrorException random_particles(medium, particle_shape, region_shape, 0.5; seed=3)

    end


end

@testset "Capsule Particle" begin

    circle_in = Circle((0.0,0.0),1.0)
    circle_out = Circle((0.0,0.0),2.0)
    a_out = Acoustic(3.0,3.0,2)
    a = Acoustic(1.0,1.0,2)
    concen_particles = [ Particle(a_out,circle_out),Particle(a,circle_in)]
    @test typeof(CapsuleParticle(concen_particles...)) <: AbstractParticle{Float64,2}

end
