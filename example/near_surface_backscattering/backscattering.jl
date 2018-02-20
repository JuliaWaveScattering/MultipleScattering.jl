using MultipleScattering
using Plots
pyplot()

listener_position = [-10.,0.]
radius = 0.5
volfrac = 0.05
max_width = 60.

bottomleft = [0.,-20.]
topright = [max_width,20.]

shape = Rectangle(bottomleft,topright)
particles = random_particles(volfrac, radius, shape; c=1.0+0.0im, ρ=0.0)

k_arr = collect(0.01:0.01:1.)
widths = 10.:10.:max_width

simulations = map(widths)  do w
    shape.topright[1] = w
    ps = filter(p -> inside(shape,p), particles)
    FrequencySimulation(ps, k_arr)
end

M = length(backscattered_waves)
differences = [norm(b - backscattered_waves[M]) for b in backscattered_waves[1:(M-1)] ]./norm(backscattered_waves[M])

plot(widths,differences)

time_arr = ω_to_t(k_arr)
time_simulations = TimeSimulation.(simulations)
plot(time_simulations...)

max_time = 40.
shape = TimeOfFlight(listener_position,max_time)

max_times = 30.:10.:max_time
M = length(max_times)
simulations = FrequencySimulation{Float64}(M)

for m in 1:M

simulations = [
    shape = TimeOfFlight(listener_position,t)
    ps = filter(p -> inside(shape,p), particles)
    FrequencySimulation(ps, k_arr)
for t in max_times]
