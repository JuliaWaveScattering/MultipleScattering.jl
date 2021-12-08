@testset "Acoustic constructors" begin
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

@testset "Acoustic special functions" begin
    T = Float64
    ms = 0:4
    map(2:3) do Dim
        L1 = basisorder_to_basislength(Acoustic{T,Dim},ms[1])
        L2 = basisorder_to_basislength(Acoustic{T,Dim},ms[end])

        m = basislength_to_basisorder(Acoustic{T,Dim}, L2)

        @test m == ms[end]
        @test L1 == 1
    end

    # Test 3D outgoing translation matrix
    ω = rand() + 0.1
    medium = Acoustic(3; ρ = 1.0, c = 1.0)
    r = rand(3) - [0.5,0.5,0.5];
    d = rand(3) - [0.5,0.5,0.5];
    d = 10 * d * norm(r) / norm(d)

    order = 2
    U = outgoing_translation_matrix(medium, 4*order, ω, d)
    vs = regular_basis_function(medium, ω)(4*order,r)
    us = outgoing_basis_function(medium, ω)(order,r + d)

    L = length(us)
    @test maximum(abs.( (U * vs)[1:L] - us) ./ abs.(us)) < 3e-6 # for linux: 4e-7

    # Test 2D outgoing translation matrix
    ω = rand() + 0.1
    medium = Acoustic(2; ρ = 1.0, c = 1.0)
    r = rand(2) - [0.5,0.5];
    d = rand(2) - [0.5,0.5];
    d = 10 * d * norm(r) / norm(d)

    # Note that to be accurate the order of vs
    order = 4
    larger_order = 3order
    U = outgoing_translation_matrix(medium, larger_order, ω, d)
    vs = regular_basis_function(medium, ω)(larger_order,r)
    us = outgoing_basis_function(medium, ω)(order,r + d)

    L = larger_order+1
    @test maximum(abs.((U * vs)[L-order:L+order] - us) ./ abs.(us)) < 1e-9
end

@testset "Acoustic circle T-matrix" begin

    circle = Sphere((0.0,0.0), 2.0)
    a2 = Acoustic(0.1, 0.1+0.0im, 2)
    a2_host = Acoustic(1.0, 1.0+0.0im, 2)
    ω = 0.5

    N_basis = estimate_outgoing_basisorder(a2_host, Particle(a2,circle), ω)

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
    p = Particle(a2, Sphere([1.0, 1.0], 0.0))
    @test_throws DomainError t_matrix(p, a2_host, ω, N_basis)

end
