# Run like a user might run it
@testset "End-to-end" begin

    @testset "Particles with same shape" begin
        circle1 = Circle((0.0,0.0),1.0)
        circle2 = Circle((0.0,5.0),2.0)
        a = Acoustic(1.0,1.0,2)
        particles = [Particle(a,circle1), Particle(a,circle2)]
        source = TwoDimAcousticPlanarSource(a,[1.0,0.0],[0.0,0.0])
        sim = FrequencySimulation(a,particles,source)
        @test true
    end

    @testset "Particles with different shape" begin
        circle = Circle((0.0,0.0),1.0)
        rect = Rectangle((2.0,2.0),3.0,2.0)
        a = Acoustic(1.0,1.0,2)
        particles = [Particle(a,circle), Particle(a,rect)]
        source = TwoDimAcousticPlanarSource(a,[1.0,0.0],[0.0,0.0])
        sim = FrequencySimulation(a,particles,source)
        @test true
    end

end
