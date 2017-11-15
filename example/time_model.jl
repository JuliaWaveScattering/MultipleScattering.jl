using MultipleScattering
using Plots

function run_time_response_single_particle(;
        k_arr = collect(linspace(0.001,1.0,1000)),
        particle_x = 100.0
    )

    # Vector of one particle
    particles = Vector{Particle{Float64}}(1)
    # Radius of particle
    radius = 1.0
    # Phase speed inside the particle
    c = 1.0 + 0.0*im
    # Density of particle
    ρ = 10.0
    # Define positions and radii
    particles[1] = Particle{Float64}([particle_x,0.0],radius,c,ρ)

    # Simulate a single particle in frequency space
    freq_model = FrequencyModel(particles,k_arr)

    # Convert the frequency model into a time model
    time_model = TimeModel(freq_model)

    return freq_model, time_model
end

function plot_time_response_single_particle()
    freq_model, time_model = run_time_response_single_particle()

    plot(
        plot(freq_model, particles),
        plot(time_model)
    )
end

"""
Test known Fourier transforms and frequency to time properties
"""
function test_fourier()
    w_arr = collect(linspace(0.0,100.0,200))
    t_arr = wTot(w_arr)
    N = length(w_arr) - 1
    t0 = t_arr[N]

    # test invertablity
    w_arr ≈ tTow(wTot(w_arr))
    fs = rand(-1:0.01:1.,2N+1,1)
    fs ≈ frequency_to_time(time_to_frequency(fs,w_arr),w_arr)

    # Gaussian to Gaussian
    as = [0.5];
    Fs = Complex{Float64}[exp(-a*w^2)*exp(im*t0*w) for w in w_arr, a in as]
    fs =  frequency_to_time(Fs,w_arr,t_arr);
    true_fs = (1/(2pi))*[sqrt(pi/a)*exp(-(t-t0)^2/(4*a)) for t in t_arr, a in as]
    fs ≈ true_fs
    Fs ≈ time_to_frequency(fs,w_arr,t_arr)
    plot(t_arr,[fs,true_fs])

    # sinc to rectangle
    rect(t) = (abs(t)<0.5) ? 1.0:0.
    as = [0.7]; #only positive values
    Fs = Complex{Float64}[ (w == 0.) ? 1:sin(a*pi*w)*exp(im*t0*w)/(a*pi*w) for w in w_arr, a in as]
    fs =  frequency_to_time(Fs,w_arr,t_arr) # only correct for half of time_arr
    true_fs = (1/(2pi))*[(1/abs(a))*rect((t-t0)/(2*pi*a)) for t in t_arr, a in as]
    mean(abs.(fs-true_fs)) < 0.02*mean(abs.(true_fs))

    # trapezoidal integration
    fs = cos.(t_arr).*exp(-(t_arr - t_arr[end]/2).^2/(t_arr[end]^2/25))
    fs = reshape(fs,length(fs),1)

    Fs_trap = time_to_frequency(fs,w_arr; method=:trapezoidal)
    Fs =time_to_frequency(fs,w_arr)
    plot(w_arr,[real(Fs),real(Fs_trap)])
    mean(abs.(Fs-Fs_trap)) < 0.0012*mean(abs.(Fs))

    fs_trap = frequency_to_time(Fs_trap,w_arr; method=:trapezoidal)
    mean(abs.(fs-fs_trap))< 1e-4*mean(abs.(fs))

    # rectangle to sinc
    # as = [1/(2*(maximum(w_arr)+w_arr[2]))]; #only positive values
    # # t0 = t_arr[end]/3
    # Fs = Complex{Float64}[ rect(a*w)*exp(im*t0*w) for w in w_arr, a in as]
    # fs =  frequency_to_time(Fs,w_arr,t_arr)
    # true_fs = [ (t==t0) ? 1/(2pi*a) : sin((t-t0)/(2a))/(pi*(t-t0)) for t in t_arr, a in as]
    # mean(abs(fs-true_fs)) < 0.04*mean(abs.(true_fs))

    # pole to wavy heaviside
    # w_arr = collect(linspace(0.0,10.0,200)) # need fine mesh near w=0
    # dw = (w_arr[2]-w_arr[1])
    # N = length(w_arr)-1;
    # time_arr = linspace(0,2π/dw,2N+2)[1:(2N+1)]
    # as = [0.1]; #only positive values
    # Fs = Complex{Float64}[1/(a-im*w) for w in w_arr, a in as]
    # fs =  frequency_to_time(Fs,w_arr,time_arr);
    # true_fs = [exp(-a*t) for t in time_arr, a in as]


    return freq_model, time_model
end
