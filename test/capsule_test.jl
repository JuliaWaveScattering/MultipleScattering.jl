import Base.Test: @testset, @test, @test_throws
import StaticArrays: SVector

using MultipleScattering

concen_particles1 = [
    Particle(Acoustic(2.0,2.0,2),Circle((0.0,0.0),1.0)),
    Particle(Acoustic(0.1,0.2,2),Circle((0.0,0.0),2.0))
]
concen_particles2 = [
    Particle(Acoustic(3.0,3.0,2),Circle((4.0,3.0),0.2)),
    Particle(Acoustic(0.4,0.8,2),Circle((4.0,3.0),1.4))
]
ps = [CapsuleParticle(concen_particles1...), CapsuleParticle(concen_particles2...)]
ps = [CapsuleParticle(concen_particles1...)]

medium = Acoustic(0.8, 0.5 + 0.0im,2)
source = TwoDimAcousticPlanarSource(medium, SVector(0.0,0.0), SVector(1.0,0.0), 1.)
sim = FrequencySimulation(medium, ps, source)

ωs = [0.1,0.2,0.3]

pressure_results, displace_results =  boundary_data(shape(ps[1].outer), ps[1].outer.medium, medium, sim, ωs; basis_order = 8)
# Continuous pressure and displacement accross particl boundary
mean(norm.(pressure_results[1].field - pressure_results[2].field))
mean(norm.(displace_results[1].field - displace_results[2].field))


ω  = 1.4

run(sim,x,ω)

bounds = Rectangle([-3.,-2.],[3.,2.]);
simres = run(sim,bounds,[ω]; res=40)

using Plots; pyplot()
plot(simres,ω,seriestype=:contour)
plot!(sim)

ωs = 0.:0.01:2.0
simres = run(sim,bounds,ωs; res=30);

ts = 0.:0.4:50.
timres = TimeSimulationResult(simres; t_vec = ts)
#timres = run(sim,bounds,0.:0.02:2.0; ts = ts, res=20)

anim = @animate for t in ts
    plot(timres,t,seriestype=:contour, clim=(0.,1.), c=:balance)
    plot!(sim)
end
#
gif(anim,"time_capsule_balance.gif", fps = 4)
# gif(anim,"time_capsule.gif", fps = 4)
