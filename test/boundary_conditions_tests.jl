@testset "boundary conditions" begin

    @testset "Circle Particle" begin

        ωs = [0.1,0.2,0.3]

        Nh = 8
        basis_order = Nh
        medium = Acoustic(1.,1.,2)

        # Choose particles
        soft_medium = Acoustic(0.,0.0 + 0.0im,2)
        p_soft = Particle(soft_medium, Sphere([1.0,2.0], .5))

        hard_medium = Acoustic(Inf,Inf + 0.0im,2)
        p_hard = Particle(hard_medium, Sphere([-3.0,-2.0], 0.3))

        sound = Acoustic(medium.ρ, 4. + 0.0im,2)
        p1 = Particle(sound, Sphere([-10.0,0.0], .2))

        particles = [p_soft, p_hard, p1]

        # Create two point sources
        source_position = SVector(0.0,0.2)
        amplitude = 1.0
        source1 = point_source(medium, source_position, amplitude);
        source2 = plane_source(medium, SVector(0.0,0.0), SVector(1.0,0.0), amplitude);
        # source2 = point_source(medium, -source_position, amplitude)
        source = 0.5*source1 + 0.5*source2;

        sim = FrequencySimulation(particles, source)
        sim_source = FrequencySimulation(source)

        result = run(sim_source, SVector(1.0,2.0), 0.1)

        pressure_results, displace_results =  boundary_data(
            shape(particles[1]), particles[1].medium, medium, sim, ωs;
            min_basis_order = 8, basis_order = 16
        )
        pressure_source_results, displace_source_results =  boundary_data(
            shape(particles[1]), particles[1].medium, medium, sim_source, ωs;
            min_basis_order = 6,
            basis_order = 10
        )

        # Zero pressure (Dirichlet) boundary condition
        @test maximum(norm.(pressure_results[1].field - pressure_results[2].field)) < 1e-6 * mean(norm.(pressure_source_results[2].field))

        pressure_results, displace_results =  boundary_data(
            shape(particles[2]), particles[2].medium, medium, sim, ωs;
            min_basis_order = 8,
            basis_order = 16
        )
        pressure_source_results, displace_source_results =  boundary_data(
            shape(particles[2]), particles[2].medium, medium, sim_source, ωs;
            basis_order = 6,
            min_basis_order = 8
        );

        # Zero displacement (Neuman) boundary condition
        @test maximum(norm.(displace_results[1].field - displace_results[2].field)) < 5e-5 * mean(norm.(displace_source_results[2].field))

        pressure_results, displace_results =  boundary_data(shape(particles[3]), particles[3].medium, medium, sim, ωs; min_basis_order = 8, basis_order = 14, dr = 8e-6);
        pressure_source_results, displace_source_results =  boundary_data(shape(particles[3]), particles[3].medium, medium, sim_source, ωs; min_basis_order = 6, basis_order = 10, dr = 1e-6);

        # Continuous pressure and displacement accross particl boundary
        @test maximum(norm.(pressure_results[1].field - pressure_results[2].field)) < 2e-6 * mean(norm.(pressure_source_results[2].field))
        @test maximum(norm.(displace_results[1].field - displace_results[2].field)) < 6e-5 * mean(norm.(displace_source_results[1].field))

        # The source pressure should always be continuous across any interface, however the displacement is only continuous because p1.medium.ρ == medium.ρ
        @test mean(norm.(pressure_source_results[1].field - pressure_source_results[2].field)) < 4e-8 * mean(norm.(pressure_source_results[2].field))
        @test mean(norm.(displace_source_results[1].field - displace_source_results[2].field)) < 5e-7 * mean(norm.(displace_source_results[1].field))
    end

    @testset "CapsuleParticle" begin

        concen_particles1 = [
            Particle(Acoustic(2.0,2.0,2),Sphere((0.0,0.0),1.0)),
            Particle(Acoustic(1.0,1.,2),Sphere((0.0,0.0),2.0))
        ]
        concen_particles2 = [
            Particle(Acoustic(2.0,1.0,2),Sphere((6.0,3.0),1.)),
            Particle(Acoustic(0.4,0.8,2),Sphere((6.0,3.0),1.4))
        ]
        particle = Particle(Acoustic(1.5,1.5,2),Sphere((-6.0,-3.0),0.8))

        ps = [CapsuleParticle(concen_particles2...), CapsuleParticle(concen_particles1...), particle]

        medium = Acoustic(0.8, 0.5 + 0.0im,2)
        source = plane_source(medium, SVector(0.0,0.0), SVector(1.0,0.0), 1.)
        sim = FrequencySimulation(ps, source)

        ωs = [0.01,0.2,0.3,1.]
        basis_vec = [8,16,16,24]

        pressure_results, displace_results =  boundary_data(shape(ps[1].outer), ps[1].outer.medium, medium, sim, ωs; basis_order_vec = basis_vec, dr = 1e9 * eps(Float64))

        # Continuous pressure and displacement accross particl boundary
        @test maximum(norm.(pressure_results[1].field - pressure_results[2].field)) / mean(norm.(pressure_results[2].field)) < 1e-6
        @test maximum(norm.(displace_results[1].field - displace_results[2].field)) / mean(norm.(displace_results[2].field)) < 5e-6

        pressure_results, displace_results =  boundary_data(shape(ps[1].inner), ps[1].inner.medium, ps[1].outer.medium, sim, ωs;  basis_order_vec = basis_vec, dr = 1e9 * eps(Float64))
        @test maximum(norm.(pressure_results[1].field - pressure_results[2].field)) / mean(norm.(pressure_results[2].field)) < 1e-6
        @test maximum(norm.(displace_results[1].field - displace_results[2].field)) / mean(norm.(displace_results[2].field)) < 5e-6

        pressure_results, displace_results =  boundary_data(shape(ps[2].outer), ps[2].outer.medium, medium, sim, ωs; basis_order_vec = basis_vec .+ 2, dr = 1e9 * eps(Float64))
        # Continuous pressure and displacement accross particl boundary
        @test maximum(norm.(pressure_results[1].field - pressure_results[2].field)) / mean(norm.(pressure_results[2].field)) < 5e-6
        @test maximum(norm.(displace_results[1].field - displace_results[2].field)) / mean(norm.(displace_results[2].field)) < 1e-5

        pressure_results, displace_results =  boundary_data(shape(ps[3]), ps[3].medium, medium, sim, ωs; basis_order_vec = basis_vec, dr = 1e9 * eps(Float64))

        @test maximum(norm.(pressure_results[1].field - pressure_results[2].field)) / mean(norm.(pressure_results[2].field)) < 1e-6
        @test maximum(norm.(displace_results[1].field - displace_results[2].field)) / mean(norm.(displace_results[2].field)) < 5e-6
    end

    @testset "Spherical Particles" begin
        # ωs = [0.1,0.2,0.3]
        ωs = [0.3]
        medium = Acoustic(3,ρ=1.1,c=0.5)

        # Choose particles
        soft_medium = Acoustic(3; ρ = 0.2, c= 0.3 + 0.0im)
        p1 = Particle(soft_medium,Sphere([1.0,2.0,0.0], .5))

        hard_medium = Acoustic(3; ρ = 3.2, c= 6.3 + 0.0im)
        p2 = Particle(hard_medium,Sphere([-1.0,-2.0,0.0], 0.3))

        sound = Acoustic(3,ρ=medium.ρ,c=1.5)
        p3 = Particle(sound,Sphere([-4.0,0.0,0.0], .2))

        particles = [p1, p2, p3]

        os = origin.(particles)
        [norm(o-u) for o in os, u in os]

        # Create two point sources
        # source_position = SVector(0.0,0.2)
        amplitude = 1.0
        # source1 = point_source(medium, source_position, amplitude);
        source2 = plane_source(medium, SVector(0.0,0.0,0.0), SVector(0.0,0.0,1.0), amplitude);
        # source2 = point_source(medium, -source_position, amplitude)
        # source = 0.5*source1 + 0.5*source2;
        source = 1.5*source2;

        sim = FrequencySimulation(particles, source)
        sim_source = FrequencySimulation(source)

        # result = run(sim_source, SVector(1.0,2.0,1.0), ωs)

        dr = 1e-7
        map(particles) do p
            pressure_results, displace_results =  boundary_data(
                shape(p), p.medium, medium, sim, ωs;
                dr = dr, basis_order = 8 #,min_basis_order = 5
            );
            pressure_source_results, displace_source_results =  boundary_data(
                shape(p), p.medium, medium, sim_source, ωs;
                dr = dr, basis_order = 8 #,min_basis_order = 5
            );

            @test maximum(norm.(pressure_results[1].field - pressure_results[2].field) ./ norm.(pressure_results[2].field)) < 1e-7
            @test maximum(norm.(displace_results[1].field - displace_results[2].field) ./ norm.(displace_results[2].field)) < 1e-6
        end
    end
end
