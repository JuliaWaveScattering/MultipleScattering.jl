using MultipleScattering

function run_time_response_single_particle(;
        ωs = LinRange(0.0,50.0,2000),
        particle_x = 100.0
    )

    # Define the particle
    radius = 1.0
    particle_medium = Acoustic(2; ρ = 10.0, c = 0.5)

    # Define positions and radii of particle
    particles = [Particle(particle_medium, Sphere([particle_x,0.0],radius))]

    # Define source wave
    host_medium = Acoustic(1.0, 1.0, 2)
    source =  plane_source(host_medium; position = [0.0,0.0], direction = [1.0,0.0])

    # Simulate a single particle in frequency space
    freq_result = run(particles, source, [0.0,0.0], ωs; min_basis_order = 1, basis_order = 14)

    # Convert the frequency simulation into a time simulation
    time_result = frequency_to_time(freq_result; impulse = GaussianImpulse(1.0; σ = 10/maximum(ωs)^2))

    return freq_result, time_result
end

function plot_time_response_single_particle()
    freq_result, time_result = run_time_response_single_particle()

    plot(
        plot(freq_result),
        plot(time_result)
    )
end
