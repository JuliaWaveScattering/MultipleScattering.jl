using MultipleScattering
using JLD
using Plots
pyplot(linewidth=2)

radius = 0.8
volfrac = 0.10
max_width = 70.

bottomleft = [0.,-max_width]
topright = [max_width,max_width]

shape = Rectangle(bottomleft,topright)
particles = random_particles(volfrac, radius, shape; c=1.0+0.0im, ρ=0.0)

listener_position = [-10.,0.]
scatter([listener_position[1]],[listener_position[2]]);
annotate!([(listener_position[1], listener_position[2] -max_width/10., "Receiver")])
plot!.(particles);
plot!(shape)

savefig("big_box.png")

k_arr = collect(0.01:0.01:1.)
widths = 10.:5.:max_width

simulations = map(widths)  do w
    shape.topright = [w,w]
    shape.bottomleft = [0,-w]
    ps = filter(p -> inside(shape,p), particles)
    FrequencySimulation(ps, k_arr)
end
save("simulations.jld", "Array{FrequencySimulation{Float64},1}",simulations)
# simulations = load("simulations.jld")["Array{FrequencySimulation{Float64},1}"];

backscattered_waves = [s.response for s in simulations]
num_particles = [length(s.particles) for s in simulations]

M = length(backscattered_waves)
bM = backscattered_waves[M] # backscattering from largest material
differences = [norm(b - bM) for b in backscattered_waves[1:(M-1)]]./norm(bM)

plot_converge = plot(num_particles[1:(M-1)], differences, xlabel = "number of particles", ylabel ="error %", label="frequency convergence")
savefig("freq_convergence.png")

time_simulations = TimeSimulation.(simulations)
widths = [round(s.shape.topright[1]) for s in simulations]
listener_position = simulations[1].listener_positions[:]
times = 2*(widths .- listener_position[1]) # time if takes for an incident plane wave to reach the furthest particles and then return to the receiver

# [Int.(round.(linspace(1,length(num_particles)-1,5))); length(num_particles)]
plot()
for i in [1,3,6,9,12,13]
    plot!(time_simulations[i],label="$(num_particles[i]) particles"
        , xlims=(0,maximum(times)+10.), ylims=(-0.6,0.3)
        , xticks = [0.; 30.; times]
    )
end
gui()
savefig("time_response.png")

time_arr = 0.:pi:34.2
time_simulations = [TimeSimulation(s;time_arr=time_arr) for s in simulations]

backscattered_waves = [s.response for s in time_simulations]
bM = backscattered_waves[M] # backscattering from largest material
differences = [norm(b - bM) for b in backscattered_waves[1:(M-1)]]./norm(bM)
plot(plot_converge)
plot!(num_particles[1:(M-1)], differences, xlabel = "number of particles", ylabel ="error %", label="time convergence")
savefig("compare_convergence.png")

near_surface_simulations = map(times) do t
    shape = TimeOfFlight(listener_position,t) # choose a material with particles only in the near surface region
    ps = filter(p -> inside(shape,p), particles) # shave off particles
    FrequencySimulation(ps, k_arr) # calculate backscattering
end
save("near_surface_simulations.jld","Array{FrequencySimulation{Float64},1}",near_surface_simulations)
