import StaticArrays: SVector
using MultipleScattering


a2 = Acoustic(0.1,0.1 + 0.0im,2)
a2_host = Acoustic(1.0,1.0 + 0.0im,2)

# Create a finite transducer by adding point sources
n = 100
amp = 10.0/n
transducer = sum(point_source(a2_host, [-2.,y], amp) for y in LinRange(-2.,2.,n))
sim = FrequencySimulation(transducer)

using Plots; gr()

bounds = Rectangle([-1.9,-5.],[10.,5.])
plot(sim, 4., seriestype=:contour, bounds=bounds, res=60)

ω_vec = 0.0:0.02:8.01
simres = run(sim, bounds, ω_vec, res=50)
timres = frequency_to_time(simres);

plot(timres, 10., seriestype=:contour)

ts = filter(t -> t<16, timres.t)
anim = @animate for t in ts
    plot(timres,t,seriestype=:contour, clim=(0.,1.25), c=:balance)
end

gif(anim,"transducer.gif", fps = 6)
