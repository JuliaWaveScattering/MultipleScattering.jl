@testset "Constructors" begin
    # 2D Acoustic
    a2 = Acoustic(0.1, 0.1+0.0im, 2)
    @test dim(a2) == 2
    @test field_dim(a2) == 1

    # 3D Acoustic
    a3 = Acoustic(1.0, 1.0+0.0im, 3)
    @test dim(a3) == 3
    @test field_dim(a3) == 1

    # Boundary condition constuctors
    @test sound_hard(2) == hard(a2)
    @test sound_hard(2) == rigid(a2)
    @test sound_hard(2) == zero_neumann(a2)

    @test sound_soft(2) == soft(a2)
    @test sound_soft(2) == pressure_release(a2)
    @test sound_soft(2) == zero_dirichlet(a2)

    @test sound_hard(2)  != sound_soft(2)
end

@testset "Circle T-matrix" begin

    circle = Circle((0.0,0.0), 2.0)
    a2 = Acoustic(0.1, 0.1+0.0im, 2)
    a2_host = Acoustic(1.0, 1.0+0.0im, 2)
    ω = 0.5
    N_basis = 5

    # Test return type is satisfied for valid input
    t = t_matrix(Particle(a2,circle), a2_host, ω, N_basis)
    @test typeof(t) <: Diagonal{Complex{Float64}}

    p = Particle(Acoustic(Inf, 0.0im, 2),circle)
    @test_throws DomainError t_matrix(p, Acoustic(1.0, 1.0+0.0im, 2), ω, N_basis)
    p = Particle(Acoustic(1.0, 1.0+0.0im, 2),circle)
    @test_throws DomainError t_matrix(p, Acoustic(0.0, Inf*im, 2), ω, N_basis)
    @test_throws DomainError t_matrix(p, Acoustic(1.0, 0.0im, 2), ω, N_basis)
    p = Particle(Acoustic(0.0, 1.0im, 2),circle)
    @test_throws DomainError t_matrix(p, Acoustic(0.0, 1.0+0.0im, 2), ω, N_basis)
    p = Particle(a2, Circle(SVector(1.0, 1.0), 0.0))
    @test_throws DomainError t_matrix(p, a2_host, ω, N_basis)

end

@testset "acoustic sources" begin
    a2 = Acoustic(0.1,0.1 + 0.0im,2)
    a2_host = Acoustic(1.0,1.0 + 0.0im,2)

    # Create two point sources
    source_position = SVector(0.0,1.0)
    amplitude = 1.0
    s1 = point_source(a2, source_position, amplitude)
    s2 = point_source(a2, 2.0*source_position, amplitude)

    # Create new souce as a linear combination of two other sources
    s3 = 2*s1 + s2

    # Check that the field is indeed a linear conbination
    x = SVector(1.0, 1.0)
    @test s3.field(x,1.0) == 2*s1.field(x,1.0) + s2.field(x,1.0)

    # Test the bessel expansions of the source
    ω = 0.8
    basis_order = 7
    centre =  SVector(1.0,0.0)

    s3_expand = source_expand(s3, centre; basis_order = 7)

    xs = [centre + 0.1 .* [cos(τ),sin(τ)] for τ = 0.0:0.3:1.5]
    @test norm([s3.field(x,ω) - s3_expand(x,ω) for x in xs]) < 1e-7*norm([s3.field(x,ω) for x in xs])

    source = plane_source(a2_host, SVector(-10.0,0.0), SVector(1.0,0.0), 1.0)
    s_expand = source_expand(source, centre; basis_order = 4)
    @test norm([source.field(x,ω) - s_expand(x,ω) for x in xs]) < 2e-9*norm([source.field(x,ω) for x in xs])
end
