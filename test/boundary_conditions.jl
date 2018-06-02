
"runs a test for the boundary conditions of penetrable particles and returns true if passed."
function boundary_conditions_test(numberofparticles::Int=4, seed = 1 )
    srand(seed) # a non-default seed may not pass all tests
    # generate 4 particles with random material properties
    particles = [Particle([0.,0.], rand(0.1:0.1:2.0), rand(0.2:0.1:10)+0.0im, rand(0.2:0.1:10)) for i=1:numberofparticles]
    shape = Rectangle(0.1, mapreduce(p->p.r,max,particles), numberofparticles)
    random_particles!(particles, shape; seed=seed) # choose random positions
    k_arr = collect(linspace(0.01,1.0,10));
    simulation = FrequencySimulation(particles,k_arr);
    boundary_data= map(4:6) do m
      simulation.hankel_order = m
      boundary_conditions(simulation,k_arr; numberofparticles=numberofparticles)
    end
    displacement_jumps = [ d[1] for d in boundary_data]
    tractions_jumps = [ d[2] for d in boundary_data]
    # to be rigorous we expect three things:
    # first, the maximum error for the displacement < 0.01 and traction < 0.01. NOTE, traction is only calculated approximately, so can be inaccurate.
    firstpass = maximum(displacement_jumps[1]) < 0.005 && maximum(tractions_jumps[1]) < 0.01
    # second, the error should, on average, increase with frequency when the hankel_order is fixed
    secondpass = issorted(mean(displacement_jumps)[:])  && issorted(mean(tractions_jumps)[:])
    # last the error should, on average, decrease when increasing the hankel_order
    lastpass = issorted(mean.(displacement_jumps), lt=(a,b) -> a>b) && issorted(mean.(tractions_jumps), lt=(a,b) -> a>b)

    return  firstpass && secondpass && lastpass
end


function boundary_fields(sim::FrequencySimulation{T,Dim,P}, ωs::Vector{T};
        numberofparticles::Int = min(4, length(sim.particles)), dr = 10000*eps(T)) where {Dim, P<:PhysicalProperties, T<:Float64}
    particles = sim.particles[1:numberofparticles]

    # points just inside particles
    inside_points = [boundary_points(p.shape; dr = - dr - 10*eps(Float64)) for p in particles]

    # points just outside particles
    outside_points = [boundary_points(p.shape; dr = 10*eps(Float64)) for p in particles]

    in_results = [run(sim, ps, ω; basis_order = basis_order) for ps in inside_points]

    x = inside1_points[1]
    map(ω->run(sim,x,ω), ωs)
    mapreduce(ω->run(sim,x,ω), union, ωs)

    [(in2_results[i].field - in1_results[i].field)/(dr * particles.medium.ρ) for i in eachindex(particles)]

    out1_results = [run(sim, ps, ω; basis_order = basis_order) for ps in outside1_points]
    out2_results = [run(sim, ps, ω; basis_order = basis_order) for ps in outside2_points]

    (out2_result.field - out1_result.field)/(dr * medium.ρ)


end

"returns (displacement_jump, stress_jump) along the boundaries of numberofparticles with wavenumbers k_arr. NOTE the stress is calculated by numerically approximately the derivative, so can be inaccurate."
function boundary_conditions(simulation::FrequencySimulation{T,Dim,P}, ωs::Vector{T};
        numberofparticles::Int = min(4, length(simulation.particles)), dr = 10000*eps(T)) where {Dim, P<:PhysicalProperties, T<:Float64}

    simulation2 = deepcopy(simulation)
    simulation2.particles = shuffle(simulation2.particles) # randomly shuffle the order
    ps = simulation2.particles[1:numberofparticles]

    positions(p,r) = hcat(points_on_boundary(p, 2*simulation.hankel_order+1; addtoradius = r)...)
    function responses(p,r)
        simulation2.listener_positions = positions(p,r)
        generate_responses!(simulation2, k_arr)
        simulation2.response
    end

    # choose listeners just outside boundary
      outside_responses = [responses(p,10*eps(T)) for p in ps]
      # listeners further out
      outside2_responses = [responses(p,dr + 10*eps(T)) for p in ps]
      diff_outside_responses = (1/simulation.ρ)*(outside2_responses - outside_responses)/dr

    # choose listeners just inside boundary
      inside_responses = [responses(p,-10*eps(T)) for p in ps]
      inside2_responses = [responses(p,-dr-10*eps(T)) for p in ps]
      diff_inside_responses = [
          (1/ps[i].ρ)*(inside_responses[i] - inside2_responses[i])/dr
      for i in eachindex(inside_responses)]

      # mean displacement and traction difference on all boundaries
      displacement = mean(mean(abs.(outside_responses .- inside_responses)),2)
      traction = mean(mean(abs.(diff_outside_responses .- diff_inside_responses)),2)

      return (displacement, traction)
end

function radial_response_simulation{T}(simulation::FrequencySimulation{T}, p::Particle{T}, k_arr::Vector{T}; opts...)
    simulation2 = deepcopy(simulation)
    # choose listeners along a radial axes
      simulation2.listener_positions = hcat(points_on_radial_axes(p; opts...)...)
      generate_responses!(simulation2, k_arr)
      return simulation2
end

"returns points on the boundary of a particle in matrix dimensions 2 x numberofpoints"
function points_on_boundary{T}(p::Particle{T}, numberofpoints::Int = 4; addtoradius::T = zero(T))
    rngθ = (0:numberofpoints-1)*2pi/numberofpoints # which angles to check for each particle
    return [p.x + (p.r+addtoradius)*[cos(θ),sin(θ)] for θ in rngθ]
end

"returns points on along a radial axes from the centre of the particle."
function points_on_radial_axes{T}(p::Particle{T}, R=T(2)*p.r, θ=zero(T); numberofpoints::Int = 20)
    return [p.x + r*[cos(θ),sin(θ)] for r in linspace(zero(T),R,numberofpoints+1)[2:end]]
end
