# Near-surface backscattering

Near-surface backscattering is a method of accurately calculating the backscattering from an infinite halfspace.
First, let us see why it is difficult to approximate the scattering from a halfspace filled with particles. That is, let us find out how many particles are required before the backscattering converges.

## Generate a large material filled with particles.

```julia
using MultipleScattering
using Plots
pyplot(linewidth=2)

radius = 0.8
volfrac = 0.10
max_width = 60.

bottomleft = [0.,-max_width]
topright = [max_width,max_width]

shape = Rectangle(bottomleft,topright)
particles = random_particles(volfrac, radius, shape; c=1.0+0.0im, ρ=0.0)
```
We will measure the backscattering at `listener_position`:

```julia
listener_position = [-10.,0.]
scatter([listener_position[1]],[listener_position[2]]);
annotate!([(listener_position[1], listener_position[2] -max_width/10., "Receiver")])
plot!.(particles);
plot!(shape)
```
![The largest quantity of particles used](big_box.png)

## Calculate backscattering for different quantity of particles
We will shave off particles on the right of this group of particles (above), and then calculate the resulting backscattered waves.
```julia
widths = 10.:10.:max_width # choose the width of the region filled with particles
k_arr = collect(0.01:0.01:1.) # choose the wavenumbers of the incident wave

simulations = map(widths) do w # this is a for loop over the array widths
    shape.topright[1] = w # choose a material with a smaller width
    ps = filter(p -> p⊆shape, particles) # shave off particles
    FrequencySimulation(ps, k_arr) # calculate backscattering
end

backscattered_waves = [s.response for s in simulations]
num_particles = [length(s.particles) for s in simulations]

M = length(backscattered_waves)
bM = backscattered_waves[M] # backscattering from largest material
differences = [norm(b - bM) for b in backscattered_waves[1:(M-1)]]./norm(bM)

plot_converge = plot(num_particles[1:(M-1)], differences, xlabel = "number of particles", ylabel ="error %", label="frequency convergence")
```
![The convergence of the response in frequency, when increasing the number of particles](freq_convergence.png)

The graph shows the rate of convergence, that is, how much the backscattering changes when including more particles (making the material deeper). The graph has not clearly converged, so we can only conclude that more than 400 particles are needed to accurately approximate the backscattering from an infinite halfspace. 

We can accelerate this convergence by considering backscattering in time.

## Calculate backscattering in time
```julia
time_simulations = TimeSimulation.(simulations)

times = 2*(widths[1:5] .- listener_position[1]) # time if takes for an incident plane wave to reach the furthest particles and then return to the receiver

plot()
for i in [1,3,6,9,12,13]
    plot!(time_simulations[i],label="$(num_particles[i]) particles"
        , xlims=(0,maximum(times)+10.), ylims=(-0.6,0.3)
        , xticks = [0.; 30.; times]
    )
end
gui()
```
![The responses in time for different quantity of particles](time_response.png)

We see that the responses in time diverge from each other more and more as time goes by. Meaning that if we only calculate the response for a short amount of time `34`, then the convergence will be accelerated.

```julia
time_arr = 0.:pi:34.2
time_simulations = [TimeSimulation(s;time_arr=time_arr) for s in simulations]

backscattered_waves = [s.response for s in time_simulations]
bM = backscattered_waves[M] # backscattering from largest material
differences = [norm(b - bM) for b in backscattered_waves[1:(M-1)]]./norm(bM)
plot(plot_converge)
plot!(num_particles[1:(M-1)], differences, xlabel = "number of particles", ylabel ="error %", label="time convergence")
```
![Compare converges for responses in time and responses in frequency](compare_convergence.png)

The convergence of the time response, for time `0<t<34`, is much faster. In fact, less than 100 particles are needed to accurately approximate the backscattering from an infinite halfspace. The reason we don't show these as log plots is because there is a small constant error (about `0.01%`) due to the discrete Fourier transform. This error is caused by the Gibbs phenomena and by assuming the backscattering is periodic (when it is not). Both these errors are well understood and can be controlled.

## Calculate backscattering only from near-surface particles
This last step is about efficiency. We want to only include particle which contribute to the backscattering for short time intervals. To do this we created a region called `TimeOfFlight(listener,time)`, where every particle in this shape takes less than `time` for their first scattered wave (due to an incident plane wave) to return to the `listener.`  More precisely, if `listener = (lx,ly)`, then every point `(x,y)` inside this shape satisfies:
`x-lx+((x-lx)^2+(y-ly)^2)^(1/2)<time` and `x>0`.

For example, look at the largest quantity of particle we used

```julia
listener_position = [-10.,0.]
shape = TimeOfFlight(listener_position,80.0)
scatter([listener_position[1]],[listener_position[2]]);
annotate!([(listener_position[1], listener_position[2] -max_width/10., "Receiver")])
plot!.(particles);
plot!(shape)
```
![Shows the particles in the shape TimeOfFlight](time_of_flight_shape.png)

For time `0<t<80` the backscattering from these particles is the same as an infinite halfspace filled with particles. To achieve this result we need only the particles inside the shape `TimeOfFlight` (region with the red outline). The particles outside this shape were unnecessary. To see this inaction:
```julia
times = 40.:15.:80.
near_surface_simulations = map(times) do t
    shape = TimeOfFlight(listener_position,t) # choose a material with particles only in the near surface region
    ps = filter(p -> p⊆shape, particles) # shave off particles
    FrequencySimulation(ps, k_arr; shape=shape) # calculate backscattering
end

time_near_simulations = TimeSimulation.(near_surface_simulations)

plot()
for i in 1:length(times)
    plot!(time_near_simulations[i],label="time of flight $(times[i])"
        , xlims=(0,maximum(times)+10.), ylims=(-0.6,0.3)
        , xticks = [0.; times], title=""
    )
end
gui()
```
![Response from particles in the shapes TimeOfFlight](time_of_flight_response.png)

Note the incident pulse has a thickness of about `10` in time, which is why the `time of flight 40` diverges from the other curves slightly before time `40`, and likewise for the other curves.
