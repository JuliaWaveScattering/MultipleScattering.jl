using MultipleScattering
using Plots

function run_time_response_single_particle(;
        k_arr = collect(LinRange(0.001,50.0,2000)),
        particle_x = 100.0
    )

    # Radius of particle
    radius = 1.0
    # Phase speed inside the particle
    c = 0.5 + 0.0*im
    # Density of particle
    ρ = 10.0
    # Define positions and radii of particle
    particles = [Particle{Float64}([particle_x,0.0],radius,c,ρ)]
    # Simulate a single particle in frequency space
    freq_simulation = FrequencySimulation(particles,k_arr; hankel_order = 10)
    # Convert the frequency simulation into a time simulation
    time_simulation = TimeSimulation(freq_simulation; impulse = get_gaussian_freq_impulse(1.0, 3./maximum(k_arr)^2))

    return freq_simulation, time_simulation
end

function plot_time_response_single_particle()
    freq_simulation, time_simulation = run_time_response_single_particle()

    plot(
        plot(freq_simulation, particles),
        plot(time_simulation)
    )
end
