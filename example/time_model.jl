using MultipleScattering
using MultipleScattering.Plot

function run_time_response_single_particle()
    k_arr = collect(linspace(0.01,1.0,100))

    # Vector of one particle
    particles = Vector{Particle{Float64}}(1)
    # Radius of particle
    radius = 1.0
    # Phase speed inside the particle
    c = 1.0 + 0.0*im
    # Density of particle
    ρ = 10.0
    # Define positions and radii
    particles[1] = Particle{Float64}([0.0,0.0],radius,c,ρ)

    # Simulate a single particle in frequency space
    freq_model = FrequencyModel(particles,k_arr)

    # Convert the frequency model into a time model
    time_model = TimeModel(freq_model)
    
    return freq_model, time_model
end

function plot_time_response_single_particle()
    freq_model, time_model = run_time_response_single_particle()

    Plots.pyplot()
    Plots.plot(plot_model(freq_model), plot_model(time_model))
end
