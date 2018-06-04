# Run like a user might run it
@testset "End-to-end" begin

    @testset "Particles with same shape" begin
        circle1 = Circle((0.0,0.0),1.0)
        circle2 = Circle((0.0,5.0),2.0)
        a_p = Acoustic(1.0,1.0,2)
        a = Acoustic(0.3,0.2,2)
        particles = [Particle(a_p,circle1), Particle(a_p,circle2)]
        source = TwoDimAcousticPlanarSource(a,[1.0,0.0],[0.0,0.0])
        sim = FrequencySimulation(a,particles,source)
        run(sim, SVector(1.0,2.0), 0.1)
        @test true
    end

    @testset "Particles with different shape" begin
        circle = Circle((0.0,0.0),1.0)
        rect = Rectangle((2.0,2.0),3.0,2.0)
        a_p = Acoustic(1.0,1.0,2)
        a = Acoustic(0.3,0.2,2)
        particles = [Particle(a_p,circle), Particle(a_p,rect)]
        source = TwoDimAcousticPlanarSource(a,[1.0,0.0],[0.0,0.0])
        sim = FrequencySimulation(a,particles,source)
        @test true
    end

end
