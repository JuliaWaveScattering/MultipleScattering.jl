import Base.Test: @testset, @test, @test_throws
import StaticArrays: SVector

using MultipleScattering

sound_p = Acoustic(.1, 0.1 + 0.0im,2)

particles = [Particle(sound_p,Circle([10.5,0.0], .5))]

sound_sim = Acoustic(1., 1. + 0.0im,2)
source = TwoDimAcousticPlanarSource(sound_sim, SVector(0.0,0.0), SVector(1.0,0.0), 1.)
sim = FrequencySimulation(sound_sim, particles, source)

ω_vec = 0.01:0.01:5.01
@test ω_vec == t_to_ω(ω_to_t(ω_vec))[2:end] # only exact for length(ω_vec) = odd number

simres = run(sim, ω_vec, SVector(0.0,0.0))

# timres = TimeSimulationResult(simres; t_vec = 0.0:0.2:50., method=:trapezoidal);
timres = TimeSimulationResult(simres; impulse = delta_freq_impulse, method=:trapezoidal);
d1 = transpose(field(timres));

timres = TimeSimulationResult(simres; method=:dft, impulse = delta_freq_impulse);
d2 = transpose(field(timres));
[d1 d2]

simres2 = FrequencySimulationResult(timres; method=:dft, impulse = delta_freq_impulse)
norm(field(simres2)[2:end] - field(simres)[:])

timres2 = TimeSimulationResult(simres2; method=:dft, impulse = delta_freq_impulse);
norm(field(timres2) - field(timres))

simres3 = FrequencySimulationResult(timres2; method=:dft, impulse = delta_freq_impulse)
norm(field(simres3) - field(simres2))

plot([real.(transpose(field(simres))), real.(transpose(field(simres2)[2:end]))])

freq_field = time_to_frequency(transpose(field(timres)), transpose(timres.t); method = :dft)
field(simres)[:] - freq_field[2:end]

plot(ω_vec, [real.(transpose(field(simres))), real.(freq_field)])

time_response = transpose(field(timres))
t_vec = transpose(timres.t)

using Plots

plot(transpose(timres.t), real.(transpose(field(timres))))
# # ω_vec
# #
field_mat = transpose(field(simres))
t_vec = ω_to_t(ω_vec)

impulse = delta_freq_impulse
impulse_vec = impulse.(transpose(ω_vec))
addzerofrequency=true
method = :dft

Dim =2; FieldDim =1; T = Float64
