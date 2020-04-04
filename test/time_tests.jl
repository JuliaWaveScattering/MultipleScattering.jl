import StaticArrays: SVector

@testset "Impulse operations" begin

    ω_vec = 0.0:0.1:1.01
    t_vec = 0.0:0.1:1.01

    gauss = GaussianImpulse(maximum(ω_vec))
    dirac = FreqDiracImpulse(ω_vec[1])
    impulse = gauss +  dirac * 2.0

    @test_throws(DomainError,[1]*dirac)

    @test all(impulse.in_freq.(ω_vec) .== gauss.in_freq.(ω_vec) + 2.0 .* dirac.in_freq.(ω_vec))
    @test all(impulse.in_time.(t_vec) .== gauss.in_time.(t_vec) + 2.0 .* dirac.in_time.(t_vec))

end

@testset "Time Result" begin
    sound_p = Acoustic(.1, 0.1 + 0.0im,2)
    particles = [Particle(sound_p,Circle([10.5,0.0], .5))]

    sound_sim = Acoustic(1., 1. + 0.0im,2)
    source = plane_source(sound_sim, [0.0,0.0], [1.0,0.0], 1.)
    sim = FrequencySimulation(particles, source)

    ω_vec = 0.0:0.01:5.01
    @test LinRange(ω_vec) == t_to_ω(ω_to_t(ω_vec)) # only exact for length(ω_vec) = even number

    # invertability of dft
        x_vec = [ [0.0,0.0], [3.0,0.0]]
        ω_vec = 0.0:0.1:1.01
        t_vec = 0.0:0.1:1.01
        simres = run(sim, x_vec, ω_vec)
        # choose an impulse which does nothing and so can be inverted
        discrete_impulse = DiscreteTimeDiracImpulse(0.0, t_vec, ω_vec)
        timres = frequency_to_time(simres; method=:dft, discrete_impulse = discrete_impulse)
        simres2 = time_to_frequency(timres; method=:dft, discrete_impulse = discrete_impulse)

        @test norm(field(simres) - field(simres2)) / norm(field(simres)) < 1e-14
        timres2 = frequency_to_time(simres2; method=:dft, impulse = TimeDiracImpulse(0.0));
        @test norm(field(timres2) - field(timres))/norm(field(timres)) < 1e-14

    # test impulse consistency by seeing if we get the same result when going from freq to time, and when going from time to freq.
        impulse = GaussianImpulse(maximum(ω_vec))
        timres2 = frequency_to_time(simres; method=:dft, impulse = impulse)
        simres2 = time_to_frequency(timres; method=:dft, impulse = impulse)
        # invert without an impulse
        timres3 = frequency_to_time(simres2; method=:dft, discrete_impulse = discrete_impulse)
        @test norm(field(timres2) - field(timres3))/norm(field(timres2)) < 1e-14

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
        # plot(timres1.t, [field(timres1)[1,:]-field(timres2)[1,:]])
end

@testset "Time shift signal" begin

    sound_sim = Acoustic(2., 5.2 + 0.0im,2)
    # suppose we measure in time the signal
    t = 0.1:0.01:30.
    t0 = 5.
    impulse_func(t) = exp( -10000.0 .* (t - t0)^4)
    measured_response = impulse_func.(t)
    plot(t, measured_response)

    ω = t_to_ω(t)
    measured_f = time_to_frequency(measured_response, t, ω)

    # at the point
    source_x = [-3.0,2.0]
    source_direction = [1.0,1.0]
    x_measure =  source_x + 4.0*source_direction

    # then can we reverse-engineer what plane wave source (of the form below) would create this measured_response
    travel_time = norm(x_measure - source_x)/real(sound_sim.c)
    impulse_time = t .- travel_time;
    impulse_response = impulse_func.(t);

    amp0 = 1.0
    source = plane_source(sound_sim, source_x, source_direction, amp0)
    sim = FrequencySimulation(source)
    simres = run(sim, [x_measure, x_measure + source_direction], ω)
    short_time = 1.:0.00005:1.5
    timres = frequency_to_time(simres; t_vec = short_time)
    i1 = findmax(abs.(field(timres)[1,:]))[2]
    i2 = findmax(abs.(field(timres)[2,:]))[2]
    # test wave peak propagation speed
    @test (short_time[i2] - short_time[i1] - norm(source_direction/sound_sim.c)) < 1e-4*norm(source_direction/sound_sim.c)

    amps = measured_f./field(simres)[1,:]

    # predict the source impulse
    predict_impulse = frequency_to_time(amps, ω, impulse_time[300:600])
    # correct mistake from ω=0
    predict_impulse = predict_impulse .- mean(predict_impulse[1:2])

    @test norm(predict_impulse - impulse_response[300:600]) < 1e-12*norm(impulse_response[300:600])
end
