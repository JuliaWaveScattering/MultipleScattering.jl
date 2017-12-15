@testset "TimeSimulation" begin
    # Time response from a single particle
    include("../example/time_response_single_particle/time_response_single_particle.jl")
    freq_model, time_model = run_time_response_single_particle()
    # Spike at start and then two diffracted spikes later, then
    # not much else

          # incident plane wave
    @test abs(time_model.response[1] - 0.98) < 1e-2     &&
          # first reflection
          abs(time_model.response[3473] - 0.041) < 1e-2 &&
          abs(time_model.time_arr[3473] - 218) < 1e-1 &&
          # largest diffraction
          abs(time_model.response[3601] - 0.0636) < 1e-3 &&
          abs(time_model.time_arr[3601] - 226) < 1e-1

    # Take samples every 200 and compare to previously run result
    @test norm(time_model.response[1:200:4000,1] - [0.9851765, -0.000428854, -0.000539422, -0.000406505, -0.000293326, -0.000361058, -0.000460193, -0.000505064, -0.000468762, -0.000440977, -0.00045894, -0.000452734, -0.000389118, -0.000326885, -0.000315177, -0.000312953, -0.000261332, -0.000161864, 0.0636417, 0.000384649]
          ) < 1e-7 # could not find how to print/showall enough digits to use ≈

    include("../example/lens/lens.jl")
    freq_model, time_model = run_lens()
    # Need to test that the spike appears at the right place
    @test true
end
