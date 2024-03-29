# # Near-surface backscattering
#
# Near-surface backscattering is a method of accurately calculating the backscattering from an infinite halfspace. For just the code see [backscattering.jl](backscattering.jl)
# First, let us see why it is difficult to approximate the scattering from a halfspace filled with particles. That is, let us find out how many particles are required before the backscattering converges.
#
# ## Generate a large material filled with particles.
#
using MultipleScattering, LinearAlgebra
using Plots

dim = 2
host_medium = Acoustic(dim; c = 1.0, ρ = 1.0)

radius = 0.8
volfrac = 0.10
max_width = 70.

particle_medium = Acoustic(dim; ρ = 0.5, c = 0.6)
particle_shape = Circle(radius)

bottomleft = [0.,-max_width]
topright = [max_width,max_width]

shape = Box([bottomleft,topright])
particles = random_particles(particle_medium, particle_shape; region_shape = shape, volume_fraction = volfrac, seed = 2)

# We send an incoming harmonic plane wave and receive the backscattering at `x`:
x = [-10.,0.]
source =  plane_source(host_medium;
    position = x,
    direction = [1.0,0.],
    amplitude = 1.0
)

plot(particles)
scatter!([x[1]],[x[2]], lab="")
annotate!([(x[1], x[2] -max_width/10., "Receiver")])
plot!(shape, linecolor = :red)

# ## Calculate backscattering for different quantity of particles
# We will shave off particles on the right of this group of particles (above), and then calculate the backscattered waves for a range of angular frequencies `ωs`.
ωs = collect(5e-3:5e-3:1.)
t_to_ω(ωs)
widths = 10.:5.0:max_width
num_particles = zeros(length(widths))

#This part below my take a while! Uncomment to run
results = map(eachindex(widths)) do i
    bottomleft = [0.,-widths[i]]
    topright = [widths[i],widths[i]]
    shape = Box([bottomleft, topright])

    ps = filter(p -> p ⊆ shape, particles) # select particles that are inside shape
    num_particles[i] = Int(length(ps))

    simulation = FrequencySimulation(ps, source)
    run(simulation, x, ωs)
end

save("results.jld2", "$(typeof(results))",results)
save("num_particles.jld2", "$(typeof(num_particles))",num_particles)

#To load results uncomment
    # results = first(values(load("results.jld2")))
    # num_particles = first(values(load("num_particles.jld2")))

backscattered_waves = field.(results)

M = length(backscattered_waves)
bM = backscattered_waves[M] # backscattering from largest material
differences = [norm(b - bM) for b in backscattered_waves[1:(M-1)]]./norm(bM)

plot_converge = plot(num_particles[1:(M-1)], differences,
    xlabel = "number of particles", yguide ="error %",
    label="frequency convergence"
)
savefig("freq_convergence.png")

gauss_impulse = GaussianImpulse(maximum(ωs))

receiver = results[1].x[1]
times = 2*(widths .- receiver[1]) # time it takes for an incident plane wave to reach the furthest particles and then return to the receiver

time_simulations = frequency_to_time.(results; impulse = gauss_impulse)

plot()
for i in [1,3,6,9,12,13]
    plot!(time_simulations[i],label="$(num_particles[i]) particles"
        , xlims=(0,maximum(times))
        , ylims=(-0.25,0.3)
        , xticks = [0.; 20.; times]
    )
end
gui()
savefig("time_response.png")

time_vec = 0.:pi:34.2
time_results = frequency_to_time.(results; t_vec = time_vec, impulse = gauss_impulse)

backscattered_waves = field.(time_results)
bM = backscattered_waves[M] # backscattering from largest material
differences = [norm(b - bM) for b in backscattered_waves[1:(M-1)]]./norm(bM)
plot(plot_converge)
plot!(num_particles[1:(M-1)], differences, xlabel = "number of particles", yguide ="error %", label="time convergence")
savefig("compare_convergence.png")

## Using near surface backscattering
shape = TimeOfFlightPlaneWaveToPoint(receiver,70.0)
scatter([receiver[1]],[receiver[2]], label = "");
annotate!([(receiver[1], receiver[2] -max_width/10., "Receiver")])
plot!(particles);
plot!(shape, linecolor=:red)

savefig("time_of_flight_shape.png")

times = 40.:10.:70.
near_surface_simulations = map(times) do t
    shape = TimeOfFlightPlaneWaveToPoint(receiver,t) # choose a material with particles only in the near surface region
    ps = filter(p -> p ⊆ shape, particles) # select particles that are inside shape
    run(FrequencySimulation(ps, source), x, ωs) # calculate backscattering
end
save("near_surface_simulations.jld2","Array{FrequencySimulation{Float64},1}",near_surface_simulations)

time_near_simulations = frequency_to_time.(near_surface_simulations; impulse = gauss_impulse)

plot()
for i in 1:length(times)
    plot!(time_near_simulations[i],label="time of flight $(times[i])"
        , xlims=(0,maximum(times)+10.)
        , ylims=(-0.2,0.3)
        , xticks = [0.; times], title=""
        , linewidth = 2.0
    )
end
gui()

savefig("time_of_flight_response.png")
