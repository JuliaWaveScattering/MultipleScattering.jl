using MultipleScattering
using Plots
pyplot()

# define two particles to scatter the wave
# the first is centered at [1.,-2.], with radius 1.0, sound speed 2. and density 10.
p1 = Particle([1.,-4.], 1.0; c = 20.0+0.0im, ρ = 10.)
p2 = Particle([3.,3.],  3.0; c = 1.0+0.0im, ρ = 0.1)
particles = [p1,p2]

# specify the angular frequency of the incident wave
w_arr = collect(0.1:0.01:1.)
# calculate and plot the frequency response
simulation = FrequencySimulation(particles, w_arr)
plot(simulation)
# the above plot used the default reciever/listener position is [-10.0; 0.0]
simulation.listener_positions
# and incident plane in the x direction
simulation.source_direction
# to change these defaults use
simulation = FrequencySimulation(particles, w_arr;
    listener_positions = [-10.,-10.],
    source_direction=[1.,1.])

# to plot the frequency response over a region that includes the particles
w = 3.2
plot(simulation,w; res=80, resp_fnc=abs)
# the green circle in the plot is the reciever position

TimeSimulation = TimeSimulation(simulation)
plot(TimeSimulation)

# simulation = FrequencySimulation(particles, w_arr; source_direction = [0.,1.], shape = Rectangle([0.,-10.],[10.,40.]))
