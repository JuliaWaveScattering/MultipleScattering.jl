# Simple random particles example

If it isn't installed, clone it from github
```julia
try using MultipleScattering
catch
    Pkg.clone("https://github.com/jondea/MultipleScattering.jl.git")
end

using MultipleScattering
```

## Define particle properties
Define the number of particles and their properties, including the region the will be placed in `shape`
```julia
num_particles = 4
radius = 1.0

particle_medium = Acoustic(0.2, 0.1, 2)
particle_shape = Circle(radius)

max_width = 20*radius
bottomleft = [0.,-max_width]
topright = [max_width,max_width]
shape = Rectangle(bottomleft,topright)

particles = random_particles(particle_medium, particle_shape, shape, num_particles)
```

Now choose the receiver position `x`, the host medium, set plane wave as a source wave, and choose the angular frequency range `ωs`
```julia
x = [-10.,0.]
host_medium = Acoustic(1.0, 1.0, 2)
source =  plane_source(host_medium; position = x, direction = [1.0,0.])

ωs = linspace(0.01,1.0,100)

simulation = FrequencySimulation(host_medium, particles, source)
result = run(simulation, x, ωs)
```

## Plot the results
We use the `Plots` package to plot both the response at the receiver position:

```julia
using Plots; pyplot(linewidth=2.0)

    plot(result, apply=real)
    plot!(result, apply=imag)
```
![Plot of response against wavenumber](plot_results.png)

And plot the whole field inside the shape `bounds` for a specific wavenumber (ω=0.8)
```julia
# plot whole field for one frequency
    bottomleft = [-15.,-max_width]
    topright = [max_width,max_width]
    bounds = Rectangle(bottomleft,topright)

    plot(simulation,0.8; res=80, bounds=bounds)
    plot!(shape, linecolor=:red)
    plot!(simulation)
    scatter!([x[1]],[x[2]], lab="receiver")    
```
![Plot real part of acoustic field](plot_field.png)

## Things to try
- Try changing the volume fraction, particle radius and k values we evaluate
