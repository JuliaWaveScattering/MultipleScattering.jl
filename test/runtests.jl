import Base.Test: @testset, @test, @test_throws
import StaticArrays: SVector

using MultipleScattering

@testset "Tests" begin
    x = SVector(1.0, 1.0)
    x2 = SVector(5.0, 5.0)
    circle = Circle(x, 2.0)
    circle_congruent = Circle(x2, 2.0)
    rect = Rectangle(x, 2.0, 3.0)

    @test volume(circle) == π*2.0^2
    @test volume(rect) == 2.0*3.0

    # 2D Acoustic
    a2 = Acoustic(0.1,0.1 + 0.0im,2)
    @test dim(a2) == 2
    @test field_dim(a2) == 1

    # 3D Acoustic
    a3 = Acoustic(1.0,1.0 + 0.0im,3)
    @test dim(a3) == 3
    @test field_dim(a3) == 1

    # Construct three particles, with two the same
    p = Particle(a2,circle)
    p_identical = Particle(a2,circle)
    p_different = Particle(a2,rect)
    p_congruent = Particle(a2,circle_congruent)

    # Test comparison operators
    @test p == p_identical
    @test p != p_different
    @test congruent(p, p_congruent)
    @test !congruent(p, p_different)

    # Cannot combine a 2D vector and shape with 3D physics
    @test_throws MethodError Particle(a3,circle)

    # Create two point sources
    source_position = SVector(0.0,1.0)
    amplitude = 1.0
    s1 = TwoDimAcousticPointSource(a2, source_position, amplitude)
    s2 = TwoDimAcousticPointSource(a2, 2.*source_position, amplitude)

    # Create new souce as a linear combination of two other sources
    s3 = 2*s1 + s2

    # Check that the field is indeed a linear conbination
    @test s3.field(x,1.0) == 2*s1.field(x,1.0) + s2.field(x,1.0)

    a2_host = Acoustic(1.0,1.0 + 0.0im,2)

    t = t_matrix(circle, a2, a2_host, 0.5, 10)
    @test typeof(t) == Diagonal{Complex{Float64}}

    @test_throws DomainError t_matrix(circle, Acoustic(Inf, 0.0im, 2), Acoustic(1.0, 1.0+0.0im, 2), 0.5, 10)
    @test_throws DomainError t_matrix(circle, Acoustic(1.0, 1.0+0.0im, 2), Acoustic(0.0, Inf*im, 2), 0.5, 10)
    @test_throws DomainError t_matrix(circle, Acoustic(1.0, 0.0im, 2), Acoustic(1.0, 0.0im, 2), 0.5, 10)
    @test_throws DomainError t_matrix(circle, Acoustic(0.0, 1.0im, 2), Acoustic(0.0, 1.0+0.0im, 2), 0.5, 10)
    @test_throws DomainError t_matrix(Circle(x, 0.0), a2, a2_host, 0.5, 10)

    # Test the bessel expansions of the source
    ω = 0.8
    centre =  SVector(1.0,0.0)
    s3_besselj = besselj_field(s3, a2, centre; basis_order = 7)
    xs = [centre + 0.1.*[cos(τ),sin(τ)] for τ = 0.0:0.3:1.5]
    @test norm([s3.field(x,ω) - s3_besselj(x,ω) for x in xs]) < 1e-7*norm([s3.field(x,ω) for x in xs])

    source = TwoDimAcousticPlanarSource(a2_host, SVector(-10.0,0.0), SVector(1.0,0.0), 1.0)
    source_besselj = besselj_field(source, a2_host, centre)
    @test norm([source.field(x,ω) - source_besselj(x,ω) for x in xs]) < 2e-9*norm([source.field(x,ω) for x in xs])

end

include("shapetests.jl")

@testset "Types" begin

    @testset "Particle" begin
        circle1 = Circle((0.0,0.0),1.0)
        circle2 = Circle((0.0,5.0),2.0)
        a = Acoustic(1.0,1.0,2)
        homog_particles = [Particle(a,circle1), Particle(a,circle2)]
        # Check types comparisons work as user would expect
        @test typeof(homog_particles) <: Particles
        @test typeof(homog_particles) <: Particles{Float64}
        @test typeof(homog_particles) <: Particles{Float64,2}

        circle = Circle((0.0,0.0),1.0)
        rect = Rectangle((2.0,2.0),3.0,2.0)
        diff_shape_particles = [Particle(a,circle), Particle(a,rect)]

        # Check types comparisons work as user would expect
        @test typeof(diff_shape_particles) <: Particles
        @test typeof(diff_shape_particles) <: Particles{Float64}
        @test typeof(diff_shape_particles) <: Particles{Float64,2}

        a2 = Acoustic(1.0,1.0,2)
        a3 = Acoustic(1.0,1.0,3)
        sphere = Sphere((0.0,0.0,0.0),1.0)

        circular_particle = Particle(a2,circle)
        @test typeof(circular_particle) <: Particle{Float64,2}

        spherical_particle = Particle(a3,sphere)
        @test typeof(spherical_particle) <: Particle{Float64,3}

        # Dimension mismatch throws error
        @test_throws MethodError Particle(a3,circle)
        @test_throws MethodError Particle(a2,sphere)

        # This is a valid vector of valid particles, but not of type Particles
        # because the dimensions don't match
        invalid_particles = [circular_particle, spherical_particle]
        @test_throws TypeError invalid_particles::Particles
    end

end

# Run like a user might run it
@testset "End-to-end" begin

    @testset "Particles with same shape" begin
        circle1 = Circle((0.0,0.0),1.0)
        circle2 = Circle((0.0,5.0),2.0)
        a = Acoustic(1.0,1.0,2)
        particles = [Particle(a,circle1), Particle(a,circle2)]
        source = TwoDimAcousticPlanarSource(a,[1.0,0.0],[0.0,0.0])
        sim = FrequencySimulation(a,particles,source)
        @test true
    end

    @testset "Particles with different shape" begin
        circle = Circle((0.0,0.0),1.0)
        rect = Rectangle((2.0,2.0),3.0,2.0)
        a = Acoustic(1.0,1.0,2)
        particles = [Particle(a,circle), Particle(a,rect)]
        source = TwoDimAcousticPlanarSource(a,[1.0,0.0],[0.0,0.0])
        sim = FrequencySimulation(a,particles,source)
        @test true
    end

end

@testset "boundary conditions" begin

    ωs = [0.1,0.2,0.3]

    Nh = 8
    basis_order = Nh
    medium = Acoustic(1.,1.,2)

    # Choose particles
    sound_soft = Acoustic(0.,0.0 + 0.0im,2)
    p_soft = Particle(sound_soft,Circle([1.0,2.0], .5))

    sound_hard = Acoustic(Inf,Inf + 0.0im,2)
    p_hard = Particle(sound_hard,Circle([-3.0,-2.0], 0.3))

    sound = Acoustic(medium.ρ, 4. + 0.0im,2)
    p1 = Particle(sound,Circle([-10.0,0.0], .2))

    particles = [p_soft, p_hard, p1]

    # Create two point sources
    source_position = SVector(0.0,0.2)
    amplitude = 1.0
    source1 = TwoDimAcousticPointSource(medium, source_position, amplitude)
    source2 = TwoDimAcousticPlanarSource(medium, SVector(0.0,0.0), SVector(1.0,0.0), amplitude)
    # source2 = TwoDimAcousticPointSource(medium, -source_position, amplitude)
    source = 0.5*source1 + 0.5*source2

    sim = FrequencySimulation(medium, particles, source)
    sim_source = FrequencySimulation(medium, source)

    pressure_results, displace_results =  boundary_data(particles[1], sim, ωs; basis_order = 8)
    pressure_source_results, displace_source_results =  boundary_data(particles[1], sim_source, ωs; basis_order = 8)

    # Zero presure (Dirichlet) boundary condition
    @test mean(norm.(pressure_results[1].field - pressure_results[2].field)) < 4e-7 * mean(norm.(pressure_source_results[2].field))

    pressure_results, displace_results =  boundary_data(particles[2], sim, ωs; basis_order = 8)
    pressure_source_results, displace_source_results =  boundary_data(particles[2], sim_source, ωs; basis_order = 8)

    # Zero displacement (Neuman) boundary condition
    @test mean(norm.(displace_results[1].field - displace_results[2].field)) < 4e-5 * mean(norm.(displace_source_results[2].field))

    pressure_results, displace_results =  boundary_data(particles[3], sim, ωs; basis_order = 8, dr = 1e-7);
    pressure_source_results, displace_source_results =  boundary_data(particles[3], sim_source, ωs; basis_order = 8, dr = 1e-7);

    # Continuous pressure and displacement accross particl boundary
    @test mean(norm.(pressure_results[1].field - pressure_results[2].field)) < 4e-9 * mean(norm.(pressure_source_results[2].field))
    @test mean(norm.(displace_results[1].field - displace_results[2].field)) < 6e-6 * mean(norm.(displace_source_results[1].field))

    # The source pressure should always be continuous accross any interface, however the displacement is only continuous because p1.medium.ρ == medium.ρ
    @test mean(norm.(pressure_source_results[1].field - pressure_source_results[2].field)) < 4e-9 * mean(norm.(pressure_source_results[2].field))
    @test mean(norm.(displace_source_results[1].field - displace_source_results[2].field)) < 5e-7 * mean(norm.(displace_source_results[1].field))

end
