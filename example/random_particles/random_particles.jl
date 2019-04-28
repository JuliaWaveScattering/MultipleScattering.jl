# If it isn't installed, clone it from github
try using MultipleScattering
catch
    Pkg.clone("https://github.com/jondea/MultipleScattering.jl.git")
end

using MultipleScattering

# Define the fraction of the volume the particles will take up, their radius and
# which wavenumbers (k) to evaluate at
num_particles = 4
radius = 1.0

particle_medium = Acoustic(2; ρ=0.2, c=0.1) # 2D particle with density ρ = 0.2 and soundspeed c = 0.1
particle_shape = Circle(radius)

max_width = 20*radius
bottomleft = [0.,-max_width]
topright = [max_width,max_width]
shape = Rectangle(bottomleft,topright)

particles = random_particles(particle_medium, particle_shape; box_shape = shape, num_particles = num_particles)

x = [-10.,0.]
host_medium = Acoustic(2; ρ=1.0, c=1.0)
source =  plane_source(host_medium; position = x, direction = [1.0,0.])

ωs = LinRange(0.01,1.0,100)

simulation = FrequencySimulation(host_medium, particles, source)
result = run(simulation, x, ωs)

# We use the `Plots` package to plot both the response at the listener position
# and the whole field for a specific wavenumber (k=0.8)
using Plots; pyplot(linewidth = 2.0)

# plot result
    plot(result, apply=real)
    plot!(result, apply=imag)

savefig("plot_result.png")

# plot whole field for one frequence
    bottomleft = [-15.,-max_width]
    topright = [max_width,max_width]
    bounds = Rectangle(bottomleft,topright)

    plot(simulation,0.8; res=80, bounds=bounds)
    plot!(shape, linecolor=:red)
    plot!(simulation)
    scatter!([x[1]],[x[2]], lab="receiver")

savefig("plot_field.png")
# ## Things to try
# - Try changing the volume fraction, particle radius and ω values we evaluate
