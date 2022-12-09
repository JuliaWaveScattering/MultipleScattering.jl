using MultipleScattering
using Plots; pyplot()
using Random
using LinearAlgebra

function box_size_convergence(  m = 4,
                                volfrac = 0.05,
                                radius = 1.0,
                                times = [20.0,30.0,40.0,50.0,60.0,70.0],
                                ω = collect(LinRange(0.01,1.0,100))
                                )

    host_medium = Acoustic(1.0, 1.0, 2)
    particle_medium = Acoustic(1.0, 0.0, 2)
    particle_shape = Circle(radius)

    listener_position = [-10.0; 0.0]
    bigshape = TimeOfFlightPlaneWaveToPoint(listener_position, maximum(times))

    seed = MersenneTwister(1).seed
    allparticles = random_particles(particle_medium, particle_shape;
                                            region_shape = bigshape,
                                            volume_fraction = volfrac,
                                            seed = seed
                                    )

    source =  plane_source( host_medium; position = [0.0, 0.0], direction = [1.0, 0.0] )

    simulations = map(times) do t
        println("Calculating response with particles at a maximum of $t away")
        shape = TimeOfFlightPlaneWaveToPoint(listener_position, t)
        particles = filter(p -> p ⊆ shape, allparticles)
        return run(particles, source, [listener_position], ω; basis_order = m)
    end

    return map(x -> frequency_to_time(x; t_vec = LinRange(0.0,100.0,201)), simulations)
end

function plot_box_size_convergence(simulations; times = [20,30,40,50,60,70])

    plot(xguide = "Time (t)", yguide ="Response")
    for s in eachindex(simulations)
        plot!(simulations[s], label = "Box cut at t=$(times[s])")
        plot!(size=(800,600))
    end
    plot!(title="Time response from random particles from increasingly large boxes")
end

simulations = box_size_convergence()
plot_box_size_convergence(simulations)
