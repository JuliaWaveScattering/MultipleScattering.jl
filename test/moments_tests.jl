@testset "Moments" begin
    begin

        # Fake results, with mean 4.0, standard deviation 2.0 and skew 0.0
        fields = [2.0, 4.0, 6.0]
        results_freq = [FrequencySimulationResult(reshape([complex(f)],1,1),[SVector(0.0)],[0.0]) for f in fields]
        results_time = [TimeSimulationResult(reshape([f],1,1),[SVector(0.0)],[0.0]) for f in fields]

        moments_freq = calculate_moments(results_freq, 3)
        @test moments_freq[1][1,1] ≈ 4.0 &&
              moments_freq[2][1,1] ≈ 2.0 &&
              moments_freq[3][1,1] ≈ 0.0

        moments_time = calculate_moments(results_time, 3)
        @test moments_time[1][1,1] ≈ 4.0 &&
              moments_time[2][1,1] ≈ 2.0 &&
              moments_time[3][1,1] ≈ 0.0
    end
end
