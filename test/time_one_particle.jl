@testset "TimeSimulation" begin
      # Time response from a single particle
      include("../docs/src/example/time_response_single_particle/time_response_single_particle.jl")
      freq_result, time_result = run_time_response_single_particle()
      # Spike at start and then two diffracted spikes later, then
      # not much else

      # incident plane wave
      @test abs( maximum(field(time_result)) - 1.0) < 1e-2

      # first diffraction
      i1 = findmin(abs.(time_result.t .- 198.5))[2]
      m, i = findmax(field(time_result)[i1-50:i1])
      t = time_result.t[i1-50:i1][i]
      @test abs(m - 0.046) < 2e-2
      @test abs(t - 198) < 1e-1

      # second diffraction
      i2 = findmin(abs.(time_result.t .- 210.0))[2]
      m, i = findmax(field(time_result)[i1+5:i2])
      t = time_result.t[i1+5:i2][i]
      @test abs(m - 0.05) < 2e-2
      @test abs(t - 206) < 2e-1

end
