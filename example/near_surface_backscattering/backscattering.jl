using MultipleScattering
using JLD
using Plots
pyplot(linewidth=2)

host_medium = Acoustic(1.0, 1.0, 2)

radius = 0.8
volfrac = 0.10
max_width = 70.

particle_medium = Acoustic(0.2, 0.1, 2)
particle_shape = Circle(radius)

bottomleft = [0.,-max_width]
topright = [max_width,max_width]

shape = Rectangle(bottomleft,topright)
particles = random_particles(particle_medium, particle_shape; box_shape = shape, volume_fraction = volfrac)

x = [-10.,0.]
source =  plane_source(host_medium; position = x,
        direction = [1.0,0.],
        amplitude = 1.0)

plot(particles)
scatter!([x[1]],[x[2]], lab="")
annotate!([(x[1], x[2] -max_width/10., "Receiver")])
plot!(shape, linecolor = :red)

savefig("big_box.png")

ωs = collect(0.01:0.01:1.)
widths = 10.:5.:max_width
num_particles = zeros(length(widths))

results = map(eachindex(widths)) do i
    bottomleft = [0.,-widths[i]]
    topright = [widths[i],widths[i]]
    shape = Rectangle(bottomleft, topright)

    ps = filter(p -> p ⊆ shape, particles) # select particles that are inside shape
    num_particles[i] = Int(length(ps))

    simulation = FrequencySimulation(host_medium, ps, source)
    run(simulation, x, ωs)
end

save("results.jld", "$(typeof(results))",results)
save("num_particles.jld", "$(typeof(num_particles))",num_particles)

# To load results uncomment
    # results = first(values(load("results.jld")))
    # num_particles = first(values(load("num_particles.jld")))

backscattered_waves = field.(results)

M = length(backscattered_waves)
bM = backscattered_waves[M] # backscattering from largest material
differences = [norm(b - bM) for b in backscattered_waves[1:(M-1)]]./norm(bM)

plot_converge = plot(num_particles[1:(M-1)], differences,
    xlabel = "number of particles", ylabel ="error %",
    label="frequency convergence"
)
savefig("freq_convergence.png")

time_simulations = frequency_to_time.(results)
receiver = results[1].x[1]
times = 2*(widths .- receiver[1]) # time it takes for an incident plane wave to reach the furthest particles and then return to the receiver

plot()
for i in [1,3,6,9,12,13]
    plot!(time_simulations[i],label="$(num_particles[i]) particles"
        , xlims=(0,maximum(times)+10.), ylims=(-0.2,0.1)
        , xticks = [0.; 30.; times]
    )
end
gui()
savefig("time_response.png")

time_vec = 0.:pi:34.2
time_results = frequency_to_time.(results; t_vec = time_vec, impulse = GaussianImpulse(maximum(ωs)))

backscattered_waves = field.(time_results)
bM = backscattered_waves[M] # backscattering from largest material
differences = [norm(b - bM) for b in backscattered_waves[1:(M-1)]]./norm(bM)
plot(plot_converge)
plot!(num_particles[1:(M-1)], differences, xlabel = "number of particles", ylabel ="error %", label="time convergence")
savefig("compare_convergence.png")

times = 40.:15.:80.
near_surface_simulations = map(times) do t
    shape = TimeOfFlight(receiver,t) # choose a material with particles only in the near surface region
    ps = filter(p -> p ⊆ shape, particles) # select particles that are inside shape
    run(FrequencySimulation(host_medium, ps, source), x, ωs) # calculate backscattering
end
save("near_surface_simulations.jld","Array{FrequencySimulation{Float64},1}",near_surface_simulations)

time_near_simulations = frequency_to_time.(near_surface_simulations; impulse = GaussianImpulse(maximum(ωs)))

plot()
for i in 1:length(times)
    plot!(time_near_simulations[i],label="time of flight $(times[i])"
        , xlims=(0,maximum(times)+10.), ylims=(-0.6,0.3)
        , xticks = [0.; times], title=""
    )
end
gui()

savefig("time_of_flight_response.png")
