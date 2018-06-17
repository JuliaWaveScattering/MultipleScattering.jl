import Base.Test: @testset, @test, @test_throws
import StaticArrays: SVector

using MultipleScattering

@testset "Time Result" begin
    sound_p = Acoustic(.1, 0.1 + 0.0im,2)
    particles = [Particle(sound_p,Circle([10.5,0.0], .5))]

    sound_sim = Acoustic(1., 1. + 0.0im,2)
    source = TwoDimAcousticPlanarSource(sound_sim, SVector(0.0,0.0), SVector(1.0,0.0), 1.)
    sim = FrequencySimulation(sound_sim, particles, source)

    ω_vec = 0.0:0.01:5.01
    @test ω_vec == t_to_ω(ω_to_t(ω_vec)) # only exact for length(ω_vec) = even number

    # invertability of dft
        x_vec = [SVector(0.0,0.0), SVector(3.0,0.0)]
        ω_vec = 0.0:0.1:1.01
        simres = run(sim, x_vec, ω_vec)
        timres = run(sim, x_vec, ω_vec; result_in_time=true, method=:dft, discrete_impulse = DiscreteTimeDiracImpulse(0.0, ω_to_t(ω_vec), ω_vec))
        simres2 = time_to_frequency(timres; method=:dft, discrete_impulse = DiscreteTimeDiracImpulse(0.0,transpose(timres.t)))
        @test norm(field(simres) - field(simres2)) / norm(field(simres)) < 1e-14
        timres2 = frequency_to_time(simres2; method=:dft, impulse = TimeDiracImpulse(0.0));
        @test norm(field(timres2) - field(timres))/norm(field(timres)) < 1e-14

    # Compare dft and trapezoidal integration
        ω_vec = 0.0:0.1:100.01 # need a high frequency to match a delta impluse function!
        simres = run(sim, x_vec, ω_vec)
        timres1 = frequency_to_time(simres; method=:trapezoidal, impulse = TimeDiracImpulse(0.0))
        timres2 = frequency_to_time(simres; method=:dft, discrete_impulse = DiscreteTimeDiracImpulse(0.0, ω_to_t(simres.ω)))
        @test norm(field(timres1) - field(timres2))/norm(field(timres1)) < 0.02

        ω_vec = 0.0:0.0001:2.01 # need a high sampling to match a delta impluse function!
        t_vec = 0.:0.5:20.
        simres = run(sim, x_vec, ω_vec)
        timres1 = frequency_to_time(simres; t_vec = t_vec, method=:trapezoidal, impulse = GaussianImpulse(maximum(simres.ω)))
        timres2 = frequency_to_time(simres; t_vec = t_vec, method=:dft, impulse = GaussianImpulse(maximum(simres.ω)))
        @test norm(field(timres1) - field(timres2))/norm(field(timres1)) < 2e-5
        # plot(timres1.t', [field(timres1)[1,:]-field(timres2)[1,:]])
end
