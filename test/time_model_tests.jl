@testset "TimeModel" begin
    # Time response from a single particle
    include("../example/time_model.jl")
    freq_model, time_model = run_time_response_single_particle()
    # Spike at start and then two diffracted spikes later, then
    # not much else

          # incident plane wave
    @test abs(time_model.response[1] - 16) < 1e-1     &&
          # first reflection
          time_model.response[3473] > 0.7             &&
          abs(time_model.time_arr[3473] - 218) < 1e-1 &&
          # largest diffraction
          time_model.response[3601] > 0.9             &&
          abs(time_model.time_arr[3601] - 226) < 1e-1

    # Take samples every 200 and compare to previously run result
    @test time_model.response[1:200:4000,1] ≈ [15.91568766260155,-0.010114700804103129,-0.010849657057909082,-0.0094274294391117,-0.008064690305764446,-0.008065776499280565,-0.008227945663215533,-0.007811684591217631,-0.0066031417092440635,-0.005390619379328917,-0.004513676793355925,-0.003388964946642101,-0.0017444062810255325,-9.962160860650897e-5,0.0011156638442981646,0.002269075560701978,0.003948284351502541,0.007065202861450666,0.9172645743771086,0.014864642396833697]

    include("../example/lens.jl")
    freq_model, time_model = run_lens()
    # Need to test that the spike appears at the right place
    @test true
end
