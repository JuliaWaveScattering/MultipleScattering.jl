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
model = FrequencyModel(particles, w_arr)
plot(model)
# the above plot used the default reciever/listener position is [-10.0; 0.0]
model.listener_positions
# and incident plane in the x direction
model.source_direction
# to change these defaults use
model = FrequencyModel(particles, w_arr;
    listener_positions = [-10.,-10.],
    source_direction=[1.,1.],
    shape = Rectangle([0.,-10.],[10.,20.]))

# to plot the frequency response over a region that includes the particles
w = 3.2
plot(model,w; res=80, resp_fnc=abs)
# the green circle in the plot is the reciever position

timemodel = TimeModel(model)
plot(timemodel)

# model = FrequencyModel(particles, w_arr; source_direction = [0.,1.], shape = Rectangle([0.,-10.],[10.,40.]))
