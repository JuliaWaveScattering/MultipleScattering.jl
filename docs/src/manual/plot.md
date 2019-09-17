# Plotting

The plotting for this package is supplied by the package [Plots](http://docs.juliaplots.org). The options and keywords used for the package Plots can be used for the plots in this package.

Below are examples of plotting the whole field in frequency (harmonic wave) and time. The examples require the package `Plots` and, mostly, `PyPlot`.

## Field - Harmonic slit
```julia
using MultipleScattering
using Plots; pyplot(size = (800,300))

radius = 1
ω = 2.0

host_medium = Acoustic(1.0, 1.0, 2)
particle_medium = Acoustic(0.0, 0.0, 2)

# Create a wall of particles
particles = [
  Particle(particle_medium, Circle([0.,y],1.0))
for y = -40:2*radius:40.]

# Make two slits in the wall
deleteat!(particles,[18,19,23,24])

# Define region to plot
bottomleft = [-10.;-15.]
topright = [30.;15.]
region = Rectangle(bottomleft, topright)

# Calculating scattering for a plane wave
source =  plane_source(host_medium; direction = [1.0,0.0])
sim = FrequencySimulation(particles, source)
result = run(sim, region, [ω]; res=100)

plot(result,ω;
    field_apply = abs, seriestype = :contour,
    title = "Absolute value"
  );
p1 = plot!(particles, ylims = (-15.0,15.0));  
plot(result,ω;
    field_apply = real, seriestype = :contour,
    title = "Real part"
);
p2 = plot!(particles, ylims = (-15.0,15.0));
plot(p1, p2)

# savefig("slit-diffraction.png")
```
![](../assets/slit-diffraction.png)

## Movie - Harmonic from random particles

```julia
using MultipleScattering
using Plots; pyplot()

num_particles = 70
radius = 1.0
ω = 1.0

host_medium = Acoustic(1.0, 1.0, 2)
particle_medium = Acoustic(0.2, 0.3, 2)
particle_shape = Circle(radius)

max_width = 50*radius
bottomleft = [0.,-max_width]
topright = [max_width,max_width]
shape = Rectangle(bottomleft,topright)

particles = random_particles(particle_medium, particle_shape; region_shape = shape, num_particles = num_particles)

source =  plane_source(host_medium; direction = [1.0,0.5])

simulation = FrequencySimulation(particles, source)

bottomleft = [-25.,-max_width]
bounds = Rectangle(bottomleft,topright)
result = run(simulation, bounds, [ω]; res=100)

ts = LinRange(0.,2pi/ω,30)

maxc = round(10*maximum(real.(field(result))))/10
minc = round(10*minimum(real.(field(result))))/10

anim = @animate for t in ts
    plot(result,ω; seriestype = :contour, phase_time=t, clim=(minc,maxc), c=:balance)
    plot!(simulation)
    plot!(colorbar=false, title="",axis=false, xlab="",ylab="")
end
#
gif(anim,"backscatter_harmonic.gif", fps = 7)
```
![backscattering from harmonic wave](../assets/backscatter_harmonic.gif)
