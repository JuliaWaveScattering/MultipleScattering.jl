import Base.Test: @testset, @test, @test_throws
import StaticArrays: SVector

using MultipleScattering

sound_p = Acoustic(1., 4. + 0.0im,2)
particles = [Particle(sound_p,Sphere([x,0.0], .5)) for x =0.:2.:10.]

sound_sim = Acoustic(1., 0.5 + 0.0im,2)
source = plane_source(sound_sim, SVector(0.0,0.0), SVector(1.0,0.0), 1.)
sim = FrequencySimulation(particles, source)

ω = 0.5

bounds = Box([[-0.,-1.], [10.,1.]])
ωs = 0.:0.1:1.
simres = run(sim, bounds, ωs)
timres = run(sim, bounds, ωs; ts = [30.], result_in_time=true)

using Plots ; #pyplot()

plot(sim,ω)
plot!(sim)

plot(simres; x = [[1,1]])
plot(simres,0.3)
plot(simres, 0.3, seriestype=:contour)

# include("plot_field.jl")
# plot_field(simres, linetype=:contour)
# plot_field(timres, linetype=:contour)
# plot!(sim)

amplitude = 1.0
# Create new souce as a linear combination of two other sources
s1 = point_source(sound_sim, SVector(8.0,0.6), amplitude)
s2 = point_source(sound_sim, SVector(2.0,-0.6), amplitude)
source = s1 + s2

sound_sim = Acoustic(1., 0.5 + 0.0im,2)
sim = FrequencySimulation(particles, source)

ω = 1.
# pyplot()
plot(sim,ω;res=60)
plot!(sim)
