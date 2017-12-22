using MultipleScattering
using Plots; pyplot()

function box_size_convergence(m=4, volfrac = 0.05,
    radius = 1.0, times = [20.0,30.0,40.0,50.0,60.0,70.0], k_arr=collect(linspace(0.01,1.0,100)) )

    listener_position = [-10.0,0.0]
    bigshape = TimeOfFlight(listener_position,maximum(times))
    
    seed = MersenneTwister(1).seed
    allparticles = random_particles(volfrac, radius, bigshape; seed = seed)

    models = map(times) do t
        println("Calculating response with particles at a maximum of $t away")
        shape = TimeOfFlight(listener_position,t)
        particles = filter(p->inside(shape,p),allparticles)
        return FrequencySimulation(particles, k_arr; seed=seed, hankel_order=m)
    end

    return map(m->TimeSimulation(m;time_arr=linspace(0.0,100.0,201)),models)
end

function plot_box_size_convergence(models;times = [20,30,40,50,60,70])

    colors = linspace(RGB(0.6,1,0.6),RGB(0,0.4,0),length(times))

    plot(xlab="Time (t)",ylab="Response")
    for m in eachindex(models)
        plot!(models[m],color=colors[m],label="Box cut at t=$(times[m])")
    end
    plot!(title="Time response from random particles from increasingly large boxes")
end

models = box_size_convergence()
plot_box_size_convergence(models)
