import StaticArrays: SVector

@testset "Constructors" begin
    # 2D Acoustic
    a2 = Acoustic(0.1, 0.1+0.0im, 2)
    @test spatial_dimension(a2) == 2
    @test field_dimension(a2) == 1

    # 3D Acoustic
    a3 = Acoustic(1.0, 1.0+0.0im, 3)
    @test spatial_dimension(a3) == 3
    @test field_dimension(a3) == 1

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
