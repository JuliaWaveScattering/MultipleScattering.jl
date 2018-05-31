import Base.Test: @testset, @test, @test_throws
import StaticArrays: SVector

using MultipleScattering

sound_p = Acoustic(.1, 0.1 + 0.0im,2)

particles = [Particle(sound_p,Circle([10.5,0.0], .5))]

sound_sim = Acoustic(1., 1. + 0.0im,2)
source = TwoDimAcousticPlanarSource(sound_sim, SVector(0.0,0.0), SVector(1.0,0.0), 1.)
sim = FrequencySimulation(sound_sim, particles, source)

ω_vec = 0.01:0.01:5.
simres = run(sim, ω_vec, SVector(0.0,0.0))

timres = TimeSimulationResult(simres; t_vec = 0.0:0.2:50.)
# timres = TimeSimulationResult(simres; impulse = delta_freq_impulse)

using Plots

plot(collect(timres.t)[:], real.(collect(field(timres))[:]))
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
