using MultipleScattering
using Plots

## Define parameters

spatial_dimension = 2
medium = Acoustic(spatial_dimension; ρ = 1.0, c = 1.0);

# define the lattice vectors
v1 = [2.0,3.0]
v2 = [-2.0,1.0]

N = 15 # 2N is the number of points along each lattice vector

radius = 0.5
particle_medium = Acoustic(spatial_dimension; ρ = 0.3, c = 0.2);

particles = [
    Particle(particle_medium, Circle(v1 * i + v2 * j, radius))
for i = -N:N, j = -2N:2N][:]

# let's create a defect by removing particles at near the origin
particles = filter(
    p -> !(p ⊆ Circle(7*radius)) && p ⊆ Circle(100*radius)
, particles)

plot(particles)

# define a radially symmetric source that is smooth
source = regular_spherical_source(medium, [1.0]; position = [0.0,0.0])
source = point_source(medium, [0.0,0.0])

# Define region to plot
bottomleft = [-50*radius;-40*radius]
topright = [50*radius;40*radius]
region = Box([bottomleft, topright])

# choose the angular frequency
ω = 1.5

# You can skip the step of defining FrequencySimulation
result = run(particles, source, region, [ω]; exclude_region = Circle(0.2), res=200)

ts = LinRange(0.,2pi/ω,30)

maxc = floor(200*maximum(real.(field(result))))/200
minc = ceil(200*minimum(real.(field(result))))/200

t = ts[1]
anim = @animate for t in ts
    plot(result,ω; seriestype = :heatmap,
        phase_time=t, clim=(minc,maxc)
        # , ylims = (-15.0,15.0) , c=:balance
    )
    plot!(particles)
    plot!(colorbar=false, title="",axis=false, xguide ="",yguide ="")
end
gif(anim,"periodic-defect-outgoing.gif", fps = 7)
