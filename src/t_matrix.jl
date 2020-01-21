"""
    t_matrix(particle, medium, ω, order)

Returns a finite T-matrix, with size depending on `order`, for a specific `particle` within a `medium` with specific physical properties.
"""
function t_matrix(p::AbstractParticle{T,Dim}, medium::PhysicalMedium{T,Dim}, ω::T, order::Integer)::AbstractMatrix{T} where {T<:AbstractFloat,Dim}

    error("T-matrix function is not yet written for $(name(p.medium)) $(name(p.shape)) in a $(name(medium)) medium")
end

get_t_matrices(medium::PhysicalMedium, species::Vector{S}, ω::AbstractFloat, Nh::Integer) where S<:Specie = get_t_matrices(medium, [s.particle for s in species], ω, Nh)

t_matrix(s::Specie, medium::PhysicalMedium, ω::AbstractFloat, order::Integer) = t_matrix(s.particle, medium, ω, order)

"""
Returns vector of T-matrices from a vector of particles in a specific domain.
Can save computation if multiple of the same kind of particle are present in the
vector.
"""
function get_t_matrices(medium::PhysicalMedium, particles::AbstractParticles, ω::AbstractFloat, Nh::Integer)::Vector

    t_matrices = Vector{AbstractMatrix}(undef,length(particles))

    # Vector of particles unique up to congruence, and the respective T-matrices
    unique_particles = Vector{AbstractParticle}(undef,0)
    unique_t_matrices = Vector{AbstractMatrix}(undef,0)

    for p_i in eachindex(particles)
        p = particles[p_i]

        # If we have calculated this T-matrix before, just point to that one
        found = false
        for cp_i in eachindex(unique_particles)
            cp = unique_particles[cp_i]
            if iscongruent(p, cp)
                t_matrices[p_i] = unique_t_matrices[cp_i]
                found = true
                break
            end
        end

        # Congruent particle was not found, we must calculate this t-matrix
        if !found
            t_matrices[p_i] = t_matrix(p, medium, ω, Nh)
            push!(unique_particles, particles[p_i])
            push!(unique_t_matrices, t_matrices[p_i])
        end

    end

    return t_matrices
end
