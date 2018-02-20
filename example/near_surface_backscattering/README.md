# Near-surface backscattering: a method of accurately calculating the backscattering from an infinite halfspace.

First, let us see why it is difficult to approximate the scattering from a halfspace filled with particles. That is, let us find out how many particles are required before the backscattering converges.

We begin by creating large box of particles.

```julia
using MultipleScattering
using Plots
pyplot()

radius = 0.5
volfrac = 0.05
max_width = 60.

bottomleft = [0.,-20.]
topright = [max_width,20.]

shape = Rectangle(bottomleft,topright)
particles = random_particles(volfrac, radius, shape; c=1.0+0.0im, ρ=0.0)
```
We will measure the backscattering at `listener_position`:

```julia
listener_position = [-10.,0.]
plot();
plot!.(particles);
scatter!(listener_position)
```

Now we will shave off particles on the right of this group of particles, and check how this effects the backscattering.
```julia
widths = 10.:10.:max_width # choose the width of the region filled with particles
k_arr = collect(0.01:0.01:1.) # choose the wavenumbers of the incident wave

simulations = map(widths) do w # this is a for loop over the array widths
    shape.topright[1] = w # choose a smaller box
    ps = filter(p -> inside(shape,p), particles) # shave off particles
    FrequencySimulation(ps, k_arr) # calculate backscattering
end

backscattered_waves =  [s.response for s in simulations]
M = length(backscattered_waves)
differences = [norm(b - backscattered_waves[M]) for b in backscattered_waves[1:(M-1)] ]./norm(backscattered_waves[M])



```
