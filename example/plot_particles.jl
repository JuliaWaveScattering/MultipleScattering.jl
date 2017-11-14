using MultipleScattering
using Plots

# We can plot a single particle...
particle = Particle([0.0,0.0])
plot(particle)

# This should work, but currently fails
# # or vector of them
# particles = random_particles(10,1.0)
# plot(particles)
