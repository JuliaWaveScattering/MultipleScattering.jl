"""
Constructor which takes volfrac, particle radius and k_arr along with a load of
optional keyword arguments. This constructor automtically generates random
particles inside a shape (default to a square Rectangle with 4 particles) then
generates the response.
"""
function FrequencySimulation{T}(volfrac::Number, radius::T, k_arr::Vector{T};
        source_direction=[one(T), zero(T)],
        c = one(Complex{T}),
        ρ = one(T),
        listener_positions = [-10one(T), zero(T)],
        source_position = [-10one(T), zero(T)],
        num_particles = 4,
        shape = Rectangle(volfrac, radius, num_particles),
        hankel_order = 3,
        seed = Base.Random.make_seed(),
        generate_responses=true
    )

    # Get the seed from the Twister so its definitely of type Vector{UInt32}
    seed = MersenneTwister(seed).seed
    particles = random_particles(volfrac, radius, shape; seed = seed)

    if isa(listener_positions,Vector)
        listener_positions = reshape(listener_positions, 2, 1)
    end

    response = Matrix{Complex{T}}(size(k_arr, 1), size(listener_positions, 2))
    simulation = FrequencySimulation{T}(
        shape, ρ, c, volfrac, particles,
        response, hankel_order,
        k_arr, listener_positions,
        source_position, source_direction,
        seed
    )
    if generate_responses generate_responses!(simulation, k_arr) end
    return simulation
end

"""
Constructor which takes a vector of particles, k_arr and a load of optional
keyword arguments. This constructor automtically generates the response.
"""
function FrequencySimulation{T}(particles::Vector{Particle{T}}, k_arr::Vector{T};
        source_direction = [one(T), zero(T)],
        c = one(Complex{T}),
        ρ = one(T),
        listener_positions = reshape([-10one(T), zero(T)], 2, 1),
        source_position = [-10one(T), zero(T)],
        shape = Rectangle(particles),
        hankel_order = 3,
        seed = Vector{UInt32}(0),
        generate_responses=true
    )
    if isa(listener_positions, Vector)
        listener_positions = reshape(listener_positions, 2, 1)
    end
    response = Matrix{Complex{T}}(size(k_arr, 1), size(listener_positions, 2))
    simulation = FrequencySimulation{T}(
        shape, ρ, c, volume(particles)/volume(shape),
        particles,
        response, hankel_order,
        k_arr, listener_positions,
        source_position, source_direction,
        seed
    )
    if generate_responses generate_responses!(simulation, k_arr) end
    return simulation
end

"Take simulation parameters, run simulation and populate the response array."
function generate_responses!{T}(simulation::FrequencySimulation{T},k_arr::Vector{T}=simulation.k_arr)
    simulation.response = Matrix{Complex{T}}(size(k_arr, 1), size(simulation.listener_positions, 2))
    # Map each k in k_arr over a the response function

    B = (2*simulation.hankel_order+1)*length(simulation.particles)
    println("\nConstructing and solving $(length(k_arr)) linear systems of size $(B)x$(B)...")
    starttime = time()
    @showprogress 0.1 "" for i=1:length(k_arr)
        simulation.response[i,:] = response(simulation,k_arr[i])
    end
    average_dur = signif((time() - starttime)/length(k_arr),3)
    println("Average time per wavelength: $average_dur seconds")
end
