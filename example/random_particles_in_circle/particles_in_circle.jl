using MultipleScattering

k_arr = collect(linspace(0.1,1.0,10))

# You can also pick your own shape, an generate random particles inside it
# with a certain radius ands volume fraction
radius = 0.3
volfrac = 0.45
centre = [0.,0.]
big_radius = 3.0

circle = Circle(big_radius,centre)
circle_simulation = FrequencySimulation(volfrac,radius,k_arr;shape=circle)

big_particle = Particle(centre,big_radius)
big_particle_simulation = FrequencySimulation([big_particle], k_arr; hankel_order=15)

using Plots
plot(
    plot(circle_simulation,0.5; drawshape=true, drawlisteners =false),
    plot(big_particle_simulation,0.5; drawlisteners =false),
    layout = (2,1)
)

plot(circle_simulation)
plot!(big_particle_simulation,title="Compare scattered wave from one big particle, \n and a circle filled with small particles")
