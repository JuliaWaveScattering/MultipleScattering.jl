import Base.Test: @testset, @test, @test_throws
import StaticArrays: SVector

using MultipleScattering

sound_p = Acoustic(1., 4. + 0.0im,2)

particles = [Particle(sound_p,Circle([5.,0.0], .5))]

sound_sim = Acoustic(1., 0.5 + 0.0im,2)
source = TwoDimAcousticPlanarSource(sound_sim, SVector(0.0,0.0), SVector(1.0,0.0), 1.)
sim = FrequencySimulation(sound_sim, particles, source)

ω_vec = 0.01:0.01:1.
simres = run(sim, ω_vec, SVector(0.0,0.0))

timres = TimeSimulationResult(simres)

field(timres)
# ω_vec
#
# field_mat = field(simres)
# t_vec = ω_to_t(ω_vec)
#
# impulse = delta_freq_impulse
# impulse_vec = impulse.(ω_vec)
# addzerofrequency=true
# method = :dft
#
# Dim =2; FieldDim =1; T = Float64
