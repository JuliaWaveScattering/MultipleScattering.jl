@testset "Basic" begin
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

    a2_host = Acoustic(1.0,1.0 + 0.0im,2)

    t = t_matrix(circle, a2, a2_host, 0.5, 10)
    @test typeof(t) == Diagonal{Complex{Float64}}

    @test_throws DomainError t_matrix(circle, Acoustic(Inf, 0.0im, 2), Acoustic(1.0, 1.0+0.0im, 2), 0.5, 10)
    @test_throws DomainError t_matrix(circle, Acoustic(1.0, 1.0+0.0im, 2), Acoustic(0.0, Inf*im, 2), 0.5, 10)
    @test_throws DomainError t_matrix(circle, Acoustic(1.0, 0.0im, 2), Acoustic(1.0, 0.0im, 2), 0.5, 10)
    @test_throws DomainError t_matrix(circle, Acoustic(0.0, 1.0im, 2), Acoustic(0.0, 1.0+0.0im, 2), 0.5, 10)
    @test_throws DomainError t_matrix(Circle(x, 0.0), a2, a2_host, 0.5, 10)

end

@testset "Particle type" begin
    circle1 = Circle((0.0,0.0),1.0)
    circle2 = Circle((0.0,5.0),2.0)
    a = Acoustic(1.0,1.0,2)
    homog_particles = [Particle(a,circle1), Particle(a,circle2)]
    # Check types comparisons work as user would expect
    @test typeof(homog_particles) <: Vector{Pt} where Pt<:AbstractParticle
    @test typeof(homog_particles) <: Vector{Pt} where Pt<:AbstractParticle{Float64}
    @test typeof(homog_particles) <: Vector{Pt} where Pt<:AbstractParticle{Float64,2}

    # CapsuleParticle should be two concentric particles
    @test_throws ErrorException CapsuleParticle(homog_particles...)

    circle_in = Circle((0.0,0.0),1.0)
    circle_out = Circle((0.0,0.0),2.0)
    a_out = Acoustic(3.0,3.0,2)
    concen_particles = [ Particle(a_out,circle_out),Particle(a,circle_in)]
    @test typeof(CapsuleParticle(concen_particles...)) <: AbstractParticle{Float64,2}

    circle = Circle((0.0,0.0),1.0)
    rect = Rectangle((2.0,2.0),3.0,2.0)
    diff_shape_particles = [Particle(a,circle), Particle(a,rect)]

    # Check types comparisons work as user would expect
    @test typeof(diff_shape_particles) <: Vector{Pt} where Pt<:AbstractParticle
    @test typeof(diff_shape_particles) <: Vector{Pt} where Pt<:AbstractParticle{Float64}
    @test typeof(diff_shape_particles) <: Vector{Pt} where Pt<:AbstractParticle{Float64,2}

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
