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
ωs = collect(LinRange(0.1,1.0,10))
result = run(simulation, x, ωs)

big_particle = Particle(particle_medium, circle)
big_particle_simulation = FrequencySimulation(host_medium, particles, source)
big_result = run(big_particle_simulation, x, ωs)

# define a bounding box for plot
    bottomleft = [-10, -2*big_radius]
    topright = [big_radius, 2*big_radius]
    box = Rectangle(bottomleft, topright)

using Plots
height = 300
pyplot(leg=false, size=(1.4*height,height))

ω = 0.5
plot(big_particle_simulation, ω; res=20, bounds = box)
p1 = plot!(big_particle)

savefig("plot_field_big.png")

plot(simulation, ω; res=20, bounds = box)
p2 = plot!(particles, linecolor = :green)

savefig("plot_field.png")

pyplot(leg=false, size=(1.8*height,height))

ωs = collect(LinRange(0.1,1.0,10))
result = run(simulation, x, ωs)
big_result = run(big_particle_simulation, x, ωs)

plot(result)
plot!(big_result, title="Compare scattered wave from one big particle, \n and a circle filled with small particles")
