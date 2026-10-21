# Run like a user might run it
@testset "End-to-end" begin

    @testset "Particles with same shape" begin

        # 2D acoustics
        circle1 = Sphere((0.0,0.0),1.0)
        circle2 = Sphere((0.0,5.0),2.0)
        a_p = Acoustic(2; ρ = 1.0, c = 1.0)
        a = Acoustic(2; ρ = 0.3, c = 0.2)
        particles = [Particle(a_p,circle1), Particle(a_p,circle2)]

        source = plane_source(a; position = [0.0,0.0], direction = [1.0,0.0])
        sim = FrequencySimulation(particles,source)

        # show(sim);

        result = run(sim, [1.0,2.0], 0.1)
        result = run(particles, source, [1.0,2.0], 0.1)
        result = 3.2*result + result*4.0im + 0.3+4.0im # changes result.field

        @test field(result)[1] == result.field[1][1] # returns
        @test field(result,1,1) == result.field[1][1] # returns

        # run simulation over a region rather than a specif set of x points
        ω = 0.1
        region = bounding_box([p.shape for p in particles])
        result = run(particles, source, region, [ω]; res=10, only_scattered_waves = true)

        # 3D acoustics
        s1 = Sphere((1.0,-1.0,2.0),1.0)
        s2 = Sphere((0.0,3.0,0.0),2.0)
        a_p = Acoustic(3; ρ = 1.0, c = 1.0)
        a = Acoustic(3; ρ = 0.3, c = 0.2)

        particles = [Particle(a_p,s1), Particle(a_p,s2)]
        source = plane_source(a; position = [0.0,0.0,0.0], direction = [0.0,0.0,1.0])

        x = [-3.0,2.0,3.0]
        ω = 0.1
        sim = FrequencySimulation(particles,source)
        result = run(particles, source, x, ω)
        result = 3.2*result + result*4.0im + 0.3+4.0im # changes result.field
        3.2 + result

        @test field(result)[1] == result.field[1][1] # returns
        @test field(result,1,1) == result.field[1][1] # returns

    end

    @testset "Particles with different shape" begin
        circle = Sphere((0.0,0.0),1.0)
        rect = Box((2.0,2.0),(3.0,2.0))
        a_p = Acoustic(1.0,1.0,2)
        a = Acoustic(0.3,0.2,2)
        particles = [Particle(a_p,circle), Particle(a_p,rect)]
        source = plane_source(a,[0.0,0.0],[1.0,0.0])
        sim = FrequencySimulation(particles,source)
        @test true
    end

    @testset "Several positions and frequencies" begin
        host = Acoustic(1.0, 1.0+0.0im, 2)
        a_p = Acoustic(0.5, 0.6+0.0im, 2)
        particles = [Particle(a_p, Sphere([0.0,0.0], 0.5)), Particle(a_p, Sphere([2.0,0.5], 0.4))]
        source = plane_source(host; direction = [1.0,0.0])
        sim = FrequencySimulation(particles, source)

        # the same number of positions and frequencies, given in any order
        xs = [[3.0,1.0], [-2.0,1.5], [0.1,0.0]]
        ωs = [1.3, 0.4, 0.9]
        result = run(sim, xs, ωs; min_basis_order = 2, basis_order = 8)
        sorted_result = run(sim, xs, sort(ωs); min_basis_order = 2, basis_order = 8)
        @test field(result)[:,1] == field(sorted_result)[:,3]
        @test field(result)[:,1] == field(run(sim, xs, ωs[1]; basis_order = 8))[:]

        f = ComplexF64[1 2; 3 4]
        @test field(FrequencySimulationResult(f, [[0.0,0.0],[1.0,0.0]], [1.0,2.0])) == f

        # PlaneSource can also be used to calculate the field inside particles
        psource = PlaneSource(host; direction = [1.0,0.0])
        x = [0.1,0.0]
        @test field(run(particles, psource, x, 0.9))[1] ≈ field(run(sim, x, 0.9))[1]

        # 3D simulation over a region, with keywords for the simulation
        host3 = Acoustic(1.0, 1.0+0.0im, 3)
        particles3 = [Particle(Acoustic(0.5, 0.6+0.0im, 3), Sphere([0.0,0.0,0.0], 0.5))]
        sim3 = FrequencySimulation(particles3, plane_source(host3; direction = [0.0,0.0,1.0]))
        result3 = run(sim3, Box([2.0,2.0,2.0]), [1.0]; res = 3, basis_order = 2)
        @test size(field(result3)) == (16,1)
    end

    include("time_one_particle.jl")

end
