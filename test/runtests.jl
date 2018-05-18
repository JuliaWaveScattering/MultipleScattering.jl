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

    ω = 0.1
    Nh = 10
    basis_order = Nh
    sound_soft = Acoustic(0.,0.1 + 0.0im,2)

    particles = [Particle(sound_soft, circle), Particle(sound_soft, circle_congruent)]
    t_matrices = get_t_matrices(a2_host, particles, ω, Nh)
    S = scattering_matrix(a2_host, particles, t_matrices, ω, Nh)

    sim = FrequencySimulation(a2_host, particles, source)


    points = boundary_points.(particles)
    listener_positions = [SVector(1.0,1.0), SVector(0.0,0.0)]
    result = run(sim, ω, listener_positions; basis_order = basis_order)
end

@testset "boundary conditions" begin

ω = 0.1
Nh = 10
basis_order = Nh
medium = Acoustic(1.,1.,2)

# Choose particles
sound_soft = Acoustic(0.,0.0 + 0.0im,2)
p_soft = Particle(sound_soft,Circle([1.0,2.0], 2.0))

sound_hard = Acoustic(Inf,100000. + 0.0im,2)
p_hard = Particle(sound_hard,Circle([-1.0,-2.0], 0.3))

sound1 = Acoustic(2.,3. + 0.0im,2)
p1 = Particle(sound1,Circle([-1.0,0.0], 0.5))

sound2 = Acoustic(4.,1.2 + 0.0im,2)
p2 = Particle(sound2,Circle([2.0,0.0], 0.5))

# t_matrix(p_hard.shape, p_hard.medium, medium, ω, Nh)
# t_matrix(p_soft.shape, p_soft.medium, medium, ω, Nh)

# Create two point sources
source_position = SVector(0.0,0.2)
amplitude = 1.0
source1 = TwoDimAcousticPointSource(medium, source_position, amplitude)
source2 = TwoDimAcousticPointSource(medium, -source_position, amplitude)
source = source1 + 1.2*source2

source = TwoDimAcousticPlanarSource(medium, SVector(-10.0,0.0), SVector(1.0,0.0), amplitude)
# This should work:
# FrequencySimulation(medium,particles,source)

particles = [p_soft, p_hard, p1]

sim = FrequencySimulation(medium, particles, source)

dr = 10000*eps(Float64)

# points just inside particles
inside1_points = [boundary_points(p.shape; dr = - dr - 10*eps(Float64)) for p in particles]
inside1_points = SVector{2,Float64}.(vcat(inside1_points...))
inside2_points = [boundary_points(p.shape; dr = - 10*eps(Float64)) for p in particles]
inside2_points = SVector{2,Float64}.(vcat(inside2_points...))

# points just outside particles
outside1_points = [boundary_points(p.shape; dr = 10*eps(Float64)) for p in particles]
outside1_points = SVector{2,Float64}.(vcat(outside1_points...))
outside2_points = [boundary_points(p.shape; dr = dr + 10*eps(Float64)) for p in particles]
outside2_points = SVector{2,Float64}.(vcat(outside2_points...))

out1_result = run(sim, ω, outside1_points; basis_order = basis_order)
out1_result.field
out2_result = run(sim, ω, outside2_points; basis_order = basis_order)
out2_result.field

(out2_result.field - out1_result.field)/(dr * medium.ρ)

in1_result = run(sim, ω, inside1_points; basis_order = basis_order)
in1_result.field
in2_result = run(sim, ω, inside2_points; basis_order = basis_order)
in2_result.field

(in2_result.field - in1_result.field)/(dr * medium.ρ)


end
