using MultipleScattering
using Plots

function run_lens(;
        a=0.25,
        volfrac=0.2,
        listener_position=[-10.0;0.0],
        outertime=35.0,
        innertime=34.0
    )

    k_arr=collect(linspace(0.01,2.0,100))

    # Generate particles which are at most outertime away from our listener
    outershape = TimeOfFlight(listener_position,outertime)
    outerparticles = random_particles(volfrac, a, outershape)

    # Filter out particles which are less than innertime away
    innershape = TimeOfFlight(listener_position,innertime)
    particles = filter(p -> !inside(innershape,p), outerparticles)

    freqsimulation = FrequencySimulation(particles, k_arr; listener_positions=listener_position)
    time_arr = 0.:0.5:90.
    return freqsimulation, TimeSimulation(freqsimulation; time_arr=time_arr, impulse = get_gaussian_freq_impulse(1.0, 4.5/maximum(k_arr)^2))
end

function plot_lens()
    simulation, time_sim = run_lens()
    plot()
    plot!.(simulation.particles);
    p = plot!();
    xticks = [0.,20.,34.,40.0,60.,80.]
    plot(p,plot(time_sim, title="", label="", xticks=xticks))
end
