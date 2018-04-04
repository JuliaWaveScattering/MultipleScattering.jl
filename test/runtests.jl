import Base.Test: @testset, @test, @test_throws

import StaticArrays: MVector

using MultipleScattering

@testset "Tests" begin
    x = MVector(1.0,1.0)
    circle = Circle(2.0)
    rect = Rectangle(2.0,3.0)

    @test volume(circle) == π*2.0^2
    @test volume(rect) == 2.0*3.0

    # 2D Acoustic
    a2 = Acoustic(1.0,1.0 + 0.0im,2)
    @test dim(a2) == 2
    @test field_dim(a2) == 1

    # 3D Acoustic
    a3 = Acoustic(1.0,1.0 + 0.0im,3)
    @test dim(a3) == 3
    @test field_dim(a3) == 1

    # Construct three particles, with two the same
    p = Particle(x,a2,circle)
    p_identical = Particle(x,a2,circle)
    p_different = Particle(x,a2,rect)

    # Test comparison operators
    @test p == p_identical
    @test p != p_different

    # Cannot combine a 2D vector and shape with 3D physics
    @test_throws MethodError Particle(x,a3,circle)

    # Create two point sources
    source_position = MVector(0.0,1.0)
    amplitude = 1.0
    s1 = TwoDimAcousticPointSource(a2, source_position, amplitude)
    s2 = TwoDimAcousticPointSource(a2, 2.*source_position, amplitude)

    # Create new souce as a linear combination of two other sources
    s3 = 2*s1 + s2

    # Check that the field is indeed a linear conbination
    @test s3.field(x,1.0) == 2*s1.field(x,1.0) + s2.field(x,1.0)

    a2_host = Acoustic(1.0,1.0 + 0.0im,2)

    t = t_matrix(circle, a2, a2_host, 0.5, 10)
    @test typeof(t) == Diagonal{Float64}

    @test_throws DomainError t_matrix(circle, Acoustic(Inf, 0.0im, 2), Acoustic(1.0, 1.0+0.0im, 2), 0.5, 10)
    @test_throws DomainError t_matrix(circle, Acoustic(1.0, 1.0+0.0im, 2), Acoustic(0.0, Inf*im, 2), 0.5, 10)
    @test_throws DomainError t_matrix(circle, Acoustic(1.0, 0.0im, 2), Acoustic(1.0, 0.0im, 2), 0.5, 10)
    @test_throws DomainError t_matrix(circle, Acoustic(0.0, 1.0im, 2), Acoustic(0.0, 1.0+0.0im, 2), 0.5, 10)
    @test_throws DomainError t_matrix(Circle(0.0), a2, a2_host, 0.5, 10)

end
