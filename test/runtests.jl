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

    ω = 0.8
    Nh = 5
    particles = [Particle(a2, circle), Particle(a2, circle_congruent)]
    t_matrices = get_t_matrices(a2_host, particles, ω, Nh)
    S = scattering_matrix(a2_host, particles, t_matrices, ω, Nh)

    source = TwoDimAcousticPlanarSource(a2_host, SVector(-10.0,0.0), SVector(1.0,0.0), 1.0)
    sim = TwoDimAcousticFrequencySimulation{Float64}(a2_host, particles, source)

    listener_positions = [SVector(1.0,1.0), SVector(0.0,0.0)]
    run(sim, ω, listener_positions)

end
