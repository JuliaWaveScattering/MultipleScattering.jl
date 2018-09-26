using MultipleScattering

# You can also pick your own shape, an generate random particles inside it
# with a certain radius ands volume fraction
radius = 0.3
volfrac = 0.45
centre = [0.,0.]
big_radius = 3.0

particle_medium = Acoustic(2; ρ=0.0, c=0.0) # 2D particle with density ρ = 0.0 and soundspeed c = 0.0
particle_shape = Circle(radius)

circle = Circle(centre, big_radius)

particles = random_particles(particle_medium, particle_shape; box_shape = circle, volume_fraction = volfrac)

plot(particles, linecolor = :green)
plot!(circle, color=:red)

x = [-10.,0.]
host_medium = Acoustic(2; ρ=1.0, c=1.0)
source =  plane_source(host_medium; position = x, direction = [1.0,0.])
simulation = FrequencySimulation(host_medium, particles, source)

# define angular frequency range
ωs = collect(linspace(0.1,1.0,10))
result = run(simulation, x, ωs)

big_particle = Particle(particle_medium, circle)
big_particle_simulation = FrequencySimulation(host_medium, particles, source)
big_result = run(big_particle_simulation, x, ωs)

# define a bounding box for plot
    bottomleft = [-10, -2*big_radius]
    topright = [big_radius, 2*big_radius]
    box = Rectangle(bottomleft, topright)

ω = 0.5
plot(big_particle_simulation, ω; res=20, bounds = box)
p1 = plot!(big_particle)

plot(simulation, ω; res=20, bounds = box)
p2 = plot!(particles, linecolor = :green)

using Plots
plot(
    p1,
    p2,
    layout = (2,1)
)

plot(circle_simulation)
plot!(big_particle_simulation,title="Compare scattered wave from one big particle, \n and a circle filled with small particles")
