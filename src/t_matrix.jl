"""
    t_matrix(particle, medium, ω, order)

Returns a finite T-matrix, with size depending on `order`, for a specific `particle` within a `medium` with specific physical properties.
"""
function t_matrix(p::AbstractParticle{Dim}, medium::PhysicalMedium{Dim}, ω, order::Integer) where {Dim}
    # ω is not typed so that any new method of t_matrix is more specific than this one
    error("T-matrix function is not yet written for a particle of type $(typeof(p)) in a medium of type $(typeof(medium)) with ω of type $(typeof(ω))")
end

"""
    get_t_matrices(PhysicalMedium, Vector{Particles}, ω, basis_order::Integer)

Returns vector of T-matrices from a vector of particles in a specific domain.
Can save computation if multiple of the same kind of particle are present in the
vector.
"""
function get_t_matrices(medium::PhysicalMedium, particles::AbstractParticles, ω::AbstractFloat, basis_order::Integer)::Vector

    t_matrices = Vector{AbstractMatrix}(undef, length(particles))

    # Vector of particles unique up to congruence, and the respective T-matrices
    unique_particles = Vector{AbstractParticle}(undef, 0)
    unique_t_matrices = Vector{AbstractMatrix}(undef, 0)

    for (i, p) in enumerate(particles)

        # If we have calculated this T-matrix before, just point to that one
        j = findfirst(cp -> iscongruent(p, cp), unique_particles)

        if isnothing(j)
            # Congruent particle was not found, we must calculate this t-matrix
            t_matrices[i] = t_matrix(p, medium, ω, basis_order)
            push!(unique_particles, p)
            push!(unique_t_matrices, t_matrices[i])
        else
            t_matrices[i] = unique_t_matrices[j]
        end
    end

    return t_matrices
end
