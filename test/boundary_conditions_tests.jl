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

    pressure_results, displace_results =  boundary_data(shape(particles[1]), particles[1].medium, medium, sim, ωs; basis_order = 8)
    pressure_source_results, displace_source_results =  boundary_data(shape(particles[1]), particles[1].medium, medium, sim_source, ωs; basis_order = 8)

    # Zero presure (Dirichlet) boundary condition
    @test mean(norm.(pressure_results[1].field - pressure_results[2].field)) < 4e-7 * mean(norm.(pressure_source_results[2].field))

    pressure_results, displace_results =  boundary_data(shape(particles[2]), particles[2].medium, medium, sim, ωs; basis_order = 8)
    pressure_source_results, displace_source_results =  boundary_data(shape(particles[2]), particles[2].medium, medium, sim_source, ωs; basis_order = 8)

    # Zero displacement (Neuman) boundary condition
    @test mean(norm.(displace_results[1].field - displace_results[2].field)) < 4e-5 * mean(norm.(displace_source_results[2].field))

    pressure_results, displace_results =  boundary_data(shape(particles[3]), particles[3].medium, medium, sim, ωs; basis_order = 8, dr = 1e-7);
    pressure_source_results, displace_source_results =  boundary_data(shape(particles[3]), particles[3].medium, medium, sim_source, ωs; basis_order = 8, dr = 1e-7);

    # Continuous pressure and displacement accross particl boundary
    @test mean(norm.(pressure_results[1].field - pressure_results[2].field)) < 4e-9 * mean(norm.(pressure_source_results[2].field))
    @test mean(norm.(displace_results[1].field - displace_results[2].field)) < 6e-6 * mean(norm.(displace_source_results[1].field))

    # The source pressure should always be continuous accross any interface, however the displacement is only continuous because p1.medium.ρ == medium.ρ
    @test mean(norm.(pressure_source_results[1].field - pressure_source_results[2].field)) < 4e-9 * mean(norm.(pressure_source_results[2].field))
    @test mean(norm.(displace_source_results[1].field - displace_source_results[2].field)) < 5e-7 * mean(norm.(displace_source_results[1].field))

end
