import Base.Test: @testset, @test, @test_throws
import StaticArrays: SVector

using MultipleScattering

concen_particles1 = [
    Particle(Acoustic(2.0,2.0,2),Circle((0.0,0.0),2.0)),
    Particle(Acoustic(0.2,0.2,2),Circle((0.0,0.0),1.0))
]
concen_particles2 = [
    Particle(Acoustic(2.0,1.0,2),Circle((6.0,3.0),1.)),
    Particle(Acoustic(0.4,0.8,2),Circle((6.0,3.0),1.4))
]
particle = Particle(Acoustic(1.5,1.5,2),Circle((-6.0,-3.0),0.8))

ps = [CapsuleParticle(concen_particles2...), CapsuleParticle(concen_particles1...), particle]
# plot(ps)

medium = Acoustic(0.8, 0.5 + 0.0im,2)
source = TwoDimAcousticPlanarSource(medium, SVector(0.0,0.0), SVector(1.0,0.0), 1.)
sim = FrequencySimulation(medium, ps, source)

ωs = [0.01,0.2,0.3,1.]

Nh = 13

pressure_results, displace_results =  boundary_data(shape(ps[1].outer), ps[1].outer.medium, medium, sim, ωs; basis_order = Nh, dr = 1e9 * eps(Float64))

# Continuous pressure and displacement accross particl boundary
maximum(norm.(pressure_results[1].field - pressure_results[2].field)) / mean(norm.(pressure_results[2].field)) < 1e-6
maximum(norm.(displace_results[1].field - displace_results[2].field)) / mean(norm.(displace_results[2].field)) < 5e-6

pressure_results, displace_results =  boundary_data(shape(ps[1].inner), ps[1].inner.medium, ps[1].outer.medium, sim, ωs; basis_order = Nh, dr = 1e9 * eps(Float64))
maximum(norm.(pressure_results[1].field - pressure_results[2].field)) / mean(norm.(pressure_results[2].field)) < 1e-6
maximum(norm.(displace_results[1].field - displace_results[2].field)) / mean(norm.(displace_results[2].field)) < 5e-6

ω  = 1.4
ps = [CapsuleParticle(concen_particles1...)]
source = TwoDimAcousticPlanarSource(medium, SVector(-3.0,0.0), SVector(1.0,0.0), 1.)
sim = FrequencySimulation(medium, ps, source)

using Plots; pyplot()

bounds = Rectangle([-3.,-2.],[3.,2.]);
simres = run(sim,bounds,[ω]; res=40, basis_order=8)

plot(simres,ω,seriestype=:contour)
plot!(sim)

ωs = 0.:0.01:2.0
# simres = run(sim,bounds,ωs; res=30, basis_order =8);
simres = run(sim,bounds,ωs; res=15, basis_order = 4);

ts = 0.:0.4:50.
timres = TimeSimulationResult(simres; t_vec = ts)
#timres = run(sim,bounds,0.:0.02:2.0; ts = ts, res=20)

t=10.0
anim = @animate for t in ts
    # plot(timres,t,seriestype=:contour, clim=(-0.8,0.8), c=:balance)
    plot(timres,t,seriestype=:contour, clim=(0.0,0.8), c=:balance)
    plot!(sim)
end
#
gif(anim,"time_capsule_balance.gif", fps = 4)
# gif(anim,"time_capsule.gif", fps = 4)


concen = CapsuleParticle(
    Particle(Acoustic(2.0,2.0,2),Circle((0.0,0.0),2.0)),
    Particle(Acoustic(0.2,0.2,2),Circle((0.0,0.0),1.0))
)

n=12
ps = [Particle(Acoustic(0.1,0.1,2),Circle((4*cos(θ),4*sin(θ)),0.5)) for θ in linspace(0.,2pi,n+1)[1:n] ]
ps = [concen; ps]
plot(ps)

ω = 1.4
source_position = SVector(0.0,-3.0)
amplitude = 1.0
source = TwoDimAcousticPointSource(medium, source_position, amplitude)
sim = FrequencySimulation(medium, ps, source)
sim_source = FrequencySimulation(medium, source)

ω = 0.01
bounds = Rectangle([-5.,-4.5],[5.,5.]);
simres = run(sim,bounds,[ω]; res=45, basis_order=12)

plot(simres,ω,seriestype=:contour)
plot!(sim)

ωs = 0.:0.01:3.0
ωs = 0.:0.01:1.0
# simres = run(sim,bounds,ωs; res=49, basis_order = 12);
# simres = run(sim,bounds,ωs; res=15, basis_order = 10);
simres = run(sim_source,bounds,ωs; res=15, basis_order = 10);

ts = 0.:0.4:50.
timres = TimeSimulationResult(simres; t_vec = ts)
#timres = run(sim,bounds,0.:0.02:2.0; ts = ts, res=20)

t=10.0
anim = @animate for t in ts
    # plot(timres,t,seriestype=:contour, clim=(-0.8,0.8), c=:balance)
    plot(timres,t,seriestype=:contour, clim=(0.0,0.8), c=:balance)
    plot!(sim)
end

gif(anim,"capsulse_dance.gif", fps = 4)
