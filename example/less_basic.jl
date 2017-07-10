using MultipleScattering
using MultipleScattering.Plot

k_Arr = collect(linspace(0.1,1.0,10))

# You can also pick your own shape, an generate random particles inside it 
# with a certain radius ands volume fraction
radius = 0.1
volfrac = 0.2

shape = Circle([0.0,0.0],5.0)
circle_model = FrequencyModel(volfrac,radius,k_arr;shape=circle)

plot_field(circle_model,0.5)


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

plot_field(two_particle_model,0.5)