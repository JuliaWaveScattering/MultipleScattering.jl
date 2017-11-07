using MultipleScattering
using MultipleScattering.Plot
import Plots

function run_lens()
    a=1.0
    volfrac=0.2
    listener_position=[-10.0;0.0]
    outertime=30.0
    innertime=20.0
    k_arr=collect(linspace(0.01,1.0,100))

    # Generate particles which are at most outertime away from our listener
    outershape = TimeOfFlight(listener_position,outertime)
    outerparticles = random_particles(volfrac, a, outershape)

    # Filter out particles which are less than innertime away
    innershape = TimeOfFlight(listener_position,innertime)
    particles = filter(p -> !inside(innershape,p), outerparticles)

    freqmodel = FrequencyModel(particles, k_arr; listener_positions=listener_position)

    return freqmodel, TimeModel(freqmodel)
end

function plot_lens()
    freqmodel, timemodel = run_time_of_flight()
    Plots.pyplot()
    Plots.plot(
        plot_particles(freqmodel),
        plot_model(timemodel)
    )
end
