using MultipleScattering
using Plots

function run_lens(;
        a=1.0,
        volfrac=0.2,
        listener_position=[-10.0;0.0],
        outertime=30.0,
        innertime=20.0
    )

    k_arr=collect(linspace(0.01,1.0,100))

    # Generate particles which are at most outertime away from our listener
    outershape = TimeOfFlight(listener_position,outertime)
    outerparticles = random_particles(volfrac, a, outershape)

    # Filter out particles which are less than innertime away
    innershape = TimeOfFlight(listener_position,innertime)
    particles = filter(p -> !inside(innershape,p), outerparticles)

    freqsimulation = FrequencySimulation(particles, k_arr; listener_positions=listener_position)

    return freqsimulation, TimeSimulation(freqsimulation; impulse = get_gaussian_freq_impulse(1.0, 3./maximum(k_arr)^2))
end

function plot_lens()
    freqsimulation, TimeSimulation = run_lens()

    plot(
        plot(freqsimulation.particles),
        plot(TimeSimulation)
    )
end
