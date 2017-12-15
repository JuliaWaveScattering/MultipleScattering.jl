using MultipleScattering

p1 = Particle([-2.,2.])
p2 = Particle([-2.,-2.])
particles = [p1,p2]

w_arr = collect(0.1:0.01:1.)
model = FrequencySimulation(particles, w_arr)

using Plots
pyplot()
plot(model)

savefig("plot_model.png")

plot(model,0.8)
savefig("plot_field.png")
