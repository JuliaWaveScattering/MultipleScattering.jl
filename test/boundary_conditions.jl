
"runs a test for the boundary conditions of penetrable particles nad returns true if passed."
function boundary_conditions_test(numberofparticles::Int=4, seed = 1 )
    srand(seed) # a non-default seed may not pass all tests
    # generate 4 particles with random material properties
    particles = [Particle([0.,0.], rand(0.1:0.1:2.0), rand(0.2:0.1:10)+0.0im, rand(0.2:0.1:10)) for i=1:numberofparticles]
    shape = Rectangle(0.1, mapreduce(p->p.r,max,particles), numberofparticles)
    random_particles!(particles, shape; seed=seed) # choose random positions
    k_arr = collect(linspace(0.01,1.0,10));
    model = FrequencyModel(particles,k_arr);
    boundary_data= map(4:6) do m
      model.hankel_order = m
      boundary_conditions(model,k_arr; numberofparticles=numberofparticles)
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


"returns (displacement_jump, stress_jump) along the boundaries of numberofparticles with wavenumbers k_arr. NOTE the stress is calculated by numerically approximately the derivative, so can be inaccurate."
function boundary_conditions{T}(model::FrequencyModel{T}, k_arr::Vector{T};
        numberofparticles::Int = min(4, length(model.particles)), dr = 100000*eps(T))

    model2 = deepcopy(model)
    model2.particles = shuffle(model2.particles) # randomly shuffle the order
    ps = model2.particles[1:numberofparticles]

    positions(p,r) = hcat(points_on_boundary(p, 2*model.hankel_order+1; addtoradius = r)...)
    function responses(p,r)
        model2.listener_positions = positions(p,r)
        generate_responses!(model2, k_arr)
        model2.response
    end

    # choose listeners just outside boundary
      outside_responses = [responses(p,10*eps(T)) for p in ps]
      # listeners further out
      outside2_responses = [responses(p,dr + 10*eps(T)) for p in ps]
      diff_outside_responses = (1/model.ρ)*(outside2_responses - outside_responses)/dr

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

function radial_response_model{T}(model::FrequencyModel{T}, p::Particle{T}, k_arr::Vector{T}; opts...)
    model2 = deepcopy(model)
    # choose listeners along a radial axes
      model2.listener_positions = hcat(points_on_radial_axes(p; opts...)...)
      generate_responses!(model2, k_arr)
      return model2
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
