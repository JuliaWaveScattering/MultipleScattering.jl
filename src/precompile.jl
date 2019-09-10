
# Precompile most objects in Float64, as it will be most commonly used
circle1 = Circle((1.0,2.0), 0.5)
circle2 = Circle((0.0,1.0), 0.5)

particle_medium = Acoustic(0.1, 0.2+0.5im, 2)
host_medium = Acoustic(1.0, 1.0+0.0im, 2)

particles = [Particle(particle_medium, circle1), Particle(particle_medium, circle2)]

source = plane_source(host_medium, SVector(-10.0,0.0), SVector(1.0,0.0), 1.0)
sim = FrequencySimulation(particles, source)

x = [SVector(1.0,1.0), SVector(0.0,0.0)]
ω = 0.5

run(sim, x, ω)
