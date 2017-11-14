using MultipleScattering

k_arr = collect(linspace(0.1,1.0,10))

# Or if you want a lot of control over the scattering problem, you can define
# each particle individually
particles = Vector{Particle{Float64}}(2)
# Phase speed inside the particle
c = 1.0 + 0.0*im
# Density of particle
ρ = 1.0
# Define their positions and radii
particles[1] = Particle{Float64}([0.2,0.8],0.1,c,ρ)
particles[2] = Particle{Float64}([0.2,-0.8],0.2,c,ρ)

two_particle_model = FrequencyModel(particles,k_arr)

using Plots
plot(two_particle_model,0.5)