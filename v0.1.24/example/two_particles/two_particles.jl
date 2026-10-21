using MultipleScattering
using Plots

# define two particles to scatter the wave
# the first is centered at [1.,-2.], with radius 1.0, sound speed 2. and density 10.
p1 = Particle(Acoustic(2; c = 20.0, ρ = 10.),Sphere([1.,-4.], 1.0))
p2 = Particle(Acoustic(2; c = 1.0, ρ = 0.1),Sphere([3.,3.], 3.0))
particles = [p1,p2]

# specify the angular frequency of the incident wave
w_arr = collect(0.1:0.01:1.)
source = plane_source(Acoustic(1.0, 1.0, 2));
# calculate and plot the frequency response at x
x = [[-10.0,0.0]];
simulation = run(particles, source, x, w_arr)
plot(simulation)
# the above plot used the reciever/listener position is [-10.0, 0.0] and incident plane wave in the default position, [0.0,0.0], and direction direction, [1.0,0.0]
# to change these defaults use
x = [[-10.0,-10.0]]
source = plane_source(Acoustic(1.0, 1.0, 2); direction = [1.0,1.0], position = [0.0,0.0]);
simulation = run(particles, source, x, w_arr)

# to plot the frequency response over a region that includes the particles
# Define region to plot
region = Box([[-11.;-11.], [6.;6.]])
# Define frequency, run the simulation and plot the field
w = 3.2
result = run(particles, source, region, [w]; res=80)
plot(result, w; field_apply=abs, seriestype = :contour)
