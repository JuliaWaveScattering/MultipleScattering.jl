@testset "Moments" begin
    begin
        # Test against a problem which can be easily solved
        models = Vector{FrequencyModel{Float64}}(3)
        # Fake responses, with mean 4.0, standard deviation 2.0 and skew 0.0
        responses = [2.0, 4.0, 6.0]
        particles = [Particle([0.0,0.0])]
        for i=1:3
            models[i] = FrequencyModel(particles,[1.0];generate_responses=false)
            models[i].response = reshape([responses[i]+0.0im],1,1)
        end
        moments = Moments(models)
        @test moments.moments[1][1] ≈ 4.0 &&
              moments.moments[2][1] ≈ 2.0 &&
              moments.moments[3][1] ≈ 0.0
    end
    begin
        # Test against a previously computed problem with a known seed
        include("../example/moments.jl")
        moments = moments_example()
        @test moments.moments[1][23] ≈ 0.9323251967911877  &&
              moments.moments[2][78] ≈ 0.2487558908409563  &&
              moments.moments[3][91] ≈ 0.15597075451712927 &&
              moments.moments[4][32] ≈ 0.17865717302214346
    end
end
