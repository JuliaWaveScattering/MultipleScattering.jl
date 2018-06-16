import Base.Test: @testset, @test, @test_throws
import StaticArrays: SVector

using MultipleScattering

sound_p = Acoustic(.1, 0.1 + 0.0im,2)

particles = [Particle(sound_p,Circle([10.5,0.0], .5))]

sound_sim = Acoustic(1., 1. + 0.0im,2)
source = TwoDimAcousticPlanarSource(sound_sim, SVector(0.0,0.0), SVector(1.0,0.0), 1.)
sim = FrequencySimulation(sound_sim, particles, source)

ω_vec = 0.0:0.01:5.01
@test ω_vec == t_to_ω(ω_to_t(ω_vec)) # only exact for length(ω_vec) = even number

x_vec = [SVector(0.0,0.0)]
x_vec = [SVector(0.0,0.0), SVector(3.0,0.0)]
ω_vec = 0.0:0.1:1.01
simres = run(sim, x_vec, ω_vec)
timres = run(sim, x_vec, ω_vec; result_in_time=true, impulse=TimeDeltaFunctionImpulse(0.0), method=:trapezoidal)

using Plots; pyplot()

bounds = Rectangle([-2.,2.], [3.,4.])
timres = run(sim, bounds, ω_vec; result_in_time=true)

plot(timres, linetype=:contour)

# timres = TimeSimulationResult(simres; t_vec = 0.0:0.2:50., method=:trapezoidal);
timres1 = frequency_to_time(simres; impulse=TimeDeltaFunctionImpulse(0.0), method=:trapezoidal)
d1 = transpose(field(timres1));

timres2 = frequency_to_time(simres; impulse=TimeDeltaFunctionImpulse(0.0), method=:dft)
d2 = transpose(field(timres2));
norm(d1 - d2)/norm(d1)

# plot(timres.t', [d1 d2])

simres2 = time_to_frequency(timres; method=:dft, impulse=FreqDeltaFunctionImpulse(0.0))
d1 = transpose(imag.(field(simres)));
d2 = transpose(imag.(field(simres2)));
# plot(simres.ω', [d1-d2])
# norm(field(simres2) - field(simres))/norm(field(simres))

timres2 = frequency_to_time(simres2; method=:dft, impulse=TimeDeltaFunctionImpulse(0.0));
norm(field(timres2) - field(timres))/norm(field(timres))

simres3 = time_to_frequency(timres2; method=:dft, impulse=FreqDeltaFunctionImpulse(0.0))
norm(field(simres3) - field(simres2))/norm(field(simres))

# plot([real.(transpose(field(simres))), real.(transpose(field(simres2)[2:end]))])

freq_field = time_to_frequency(transpose(field(timres)), transpose(timres.t); method = :dft)
# field(simres)[:] - freq_field[2:end]

# plot(ω_vec, [real.(transpose(field(simres))), real.(freq_field)])

time_response = transpose(field(timres))
t_vec = transpose(timres.t)

using Plots

plot(transpose(timres.t), real.(transpose(field(timres))))
# # ω_vec
# #
field_mat = transpose(field(simres))
t_vec = ω_to_t(ω_vec)

impulse = TimeDeltaFunctionImpulse(0.0)
impulse_in_freq_vec = impulse.in_freq.(transpose(ω_vec))
addzerofrequency=true
method = :dft

Dim =2; FieldDim =1; T = Float64
