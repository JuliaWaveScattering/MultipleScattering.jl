using MultipleScattering

p1 = Particle([-2.,2.])
p2 = Particle([-2.,-2.])
particles = [p1,p2]

w_arr = collect(0.1:0.01:1.)
simulation = FrequencySimulation(particles, w_arr)

using Plots
pyplot()
plot(simulation)

savefig("plot_simulation.png")

plot(simulation,0.8)
savefig("plot_field.png")
