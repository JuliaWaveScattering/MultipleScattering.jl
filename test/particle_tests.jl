@testset "Particle" begin

    a = Acoustic(1.0, 1.0, 2)
    a2 = Acoustic(1.0, 2.0, 2)

    @testset "Types" begin

        circle1 = Sphere((0.0,0.0),1.0)
        circle2 = Circle((0.0,5.0),2.0)
        a = Acoustic(1.0,1.0,2)
        homog_particles = [Particle(a,circle1), Particle(a,circle2)]

        # Check types comparisons work as user would expect
        @test typeof(homog_particles) <: AbstractParticles
        @test typeof(homog_particles) <: AbstractParticles{2}

        circle = Sphere((0.0,0.0),1.0)
        rect = Box((2.0,2.0),(3.0,2.0))
        diff_shape_particles = [Particle(a,circle), Particle(a,rect)]

        # Check types comparisons work as user would expect
        @test typeof(diff_shape_particles) <: AbstractParticles
        @test typeof(diff_shape_particles) <: AbstractParticles{2}

        a2 = Acoustic(1.0,1.0,2)
        a3 = Acoustic(1.0,1.0,3)
        sphere = Sphere((0.0,0.0,0.0),1.0)

        circular_particle = Particle(a2,circle)
        @test typeof(circular_particle) <: Particle{2}

        spherical_particle = Particle(a3,sphere)
        @test typeof(spherical_particle) <: AbstractParticle{3}

        # Dimension mismatch throws error
        @test_throws MethodError Particle(a3,circle)
        @test_throws MethodError Particle(a2,sphere)

        # This is a valid vector of valid particles, but not of type Particles
        # because the dimensions don't match
        invalid_particles = [circular_particle, spherical_particle]
        @test_throws TypeError invalid_particles::(Vector{Pt} where {Dim, Pt <: AbstractParticle{Dim}})

        # does not throw an error
        diff_shape_particles::(Vector{Pt} where {Dim, Pt <: AbstractParticle{Dim}})
        @test true

    end

    @testset "Comparisons" begin

        circle = Sphere((1.0,3.0),2.0)
        circle_identical = Sphere((1.0,3.0),2.0)
        circle_congruent = Sphere((4.0,7.0),2.0)
        rect = Box((1.0,2.0),(2.0,4.0))

        resonator = SphericalHelmholtz((1.0,3.0),2.0, 0.2, 0.01, -1.3)
        resonator_kws = SphericalHelmholtz((1.0,3.0),2.0; inner_radius = 0.2, aperture = 0.1, orientation = -1.3)
        resonator_dif_aperture = SphericalHelmholtz((1.0,3.0),2.0; aperture = 0.1)
        resonator_identical = SphericalHelmholtz((1.0,3.0), 2.0; orientation = -1.3, inner_radius = 0.2, aperture = 0.01)
        resonator_congruent = SphericalHelmholtz((4.0,7.0), 2.0; orientation = -1.3, inner_radius = 0.2, aperture = 0.01)

        # Construct four particles, with two the same
        p = Particle(a2,circle)
        p_reference = p
        p_identical = Particle(a2,circle_identical)
        p_different = Particle(a2, rect)
        p_congruent = Particle(a2,circle_congruent)

        # Construct three resonator particles
        pr = Particle(a2, resonator)
        pr_reference = pr
        pr_dif_aperture = Particle(a2, resonator_dif_aperture)
        pr_identical = Particle(a2, resonator_identical)
        pr_different = p_different
        pr_congruent = Particle(a2, resonator_congruent)

        # Test comparison operators
        @test p == p_identical
        @test p != p_different
        @test !(p == p_different)
        @test pr != pr_dif_aperture
        @test !(pr == pr_dif_aperture)
        @test p != pr
        @test !(p == pr)
        @test iscongruent(p, p_congruent)
        @test !iscongruent(p, p_different)
        @test !iscongruent(p, pr)
        @test pr == pr_identical
        @test pr != pr_different
        @test !(pr == pr_different)
        @test iscongruent(pr, pr_congruent)
        @test !iscongruent(pr, pr_dif_aperture)
        @test !iscongruent(pr, pr_different)

        # Check that Julia fallback === works
        @test p === p_reference
        @test pr === pr_reference
        # Particles are immutable, so are compared at the bit levelS
        @test p === p_identical
        @test pr === pr_identical

        # Construct two almost identical particles
        p1 = Particle(a, Sphere((0.0,0.0),1.0))
        p1b = Particle(a, Sphere((-0.0,0.0),1.0))

        # isequal is strict about +/-0.0 but IEEE compliant == is not strict
        @test p1 == p1b
        @test !isequal(p1,p1b)

    end

    @testset "Plot" begin
        # Just run it to see if we have any errors (yes thats a very low bar)

        # Single particle
        particle = Particle(a, Sphere((1.0,0.0), 2.0))
        plot(particle)

        # Vector of particles of same shape
        particles = [
            Particle(a, Sphere((1.0,0.0), 2.0)),
            Particle(a, Sphere((3.0,-1.0), 1.0))
        ]
        plot(particles)

        # Vector of particles of different shape
        particles = [
            CapsuleParticle(
                Particle(Acoustic(2.0,2.0,2),Sphere((5.0,3.0),1.0)),
                Particle(Acoustic(1.0,1.,2),Sphere((5.0,3.0),2.0))
            ),
            Particle(a, Box((1.0,0.0),(2.0, 3.0))),
            Particle(a, Sphere((3.0,-1.0), 1.0))
        ]
        plot(particles)

        @test true
    end

end

@testset "Capsule Particle" begin

    circle_in = Sphere((0.0,0.0),1.0)
    circle_out = Sphere((0.0,0.0),2.0)
    a_out = Acoustic(3.0,3.0,2)
    a = Acoustic(1.0,1.0,2)
    concen_particles = [ Particle(a_out,circle_out),Particle(a,circle_in)]
    @test typeof(CapsuleParticle(concen_particles...)) <: AbstractParticle{2}

end
