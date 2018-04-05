
"Create the matrix S we invert to find the coefficients based on the forcing"
function scattering_matrix(medium::PhysicalProperties, particles::Vector, ω::T, Nh::Integer)::Matrix{Complex{T}} where T
    # Generate response for one specific k
    # Number of particles
    P = length(particles)

    # No particles means no scattering
    if P == 0
        warning("You have computed the scattering matrix with no particles, are you sure something hasn't gone wrong?")
        return Matrix{Complex{T}}(0,0)
    end

    # Number of hankel basis function at each particle
    H = 2Nh + 1

    t_matrices = get_t_matrices(medium, particles, ω, Nh)

    basis_function = get_basis_function(medium, ω)

    # Faire: this could potentially return an MMatrix
    function S_block(j,l)
        if j == l
            return eye(Complex{T}, H, H)
        else
            x_jl = particles[j].position .- particles[l].position
            # Faire: basis functions could be more efficient if it returned a vector
            basis_vec = OffsetArray(map(m->basis_function(m,x_jl), -2Nh:2Nh),-2Nh:2Nh)
            mat = [basis_vec[p-m] for m in -Nh:Nh, p in -Nh:Nh]
            return -t_matrices[j] * mat
        end
    end

    # Evaluates all the blocks in the big matrix
    S_blocks = [S_block(j,l) for j in 1:P, l in 1:P]

    # Reshape S_blocks into big matrix
    S = Matrix{Complex{T}}(P*H, P*H)
    for i in 1:P
        for j in 1:P
            S[(i-1)*H+1:i*H, (j-1)*H+1:j*H] .= S_blocks[i,j]
        end
    end

    return S
end

function get_t_matrices(medium::PhysicalProperties, particles::Vector, ω::AbstractFloat, Nh::Integer)::Vector

    t_matrices = Vector{AbstractMatrix}(length(particles))
    congruent_particles = Vector{Tuple{PhysicalProperties,Shape}}(0)

    for p_i in eachindex(particles)
        p = particles[p_i]

        # If we have calculated this T-matrix before, just link to that one
        found = false
        for cp_i in eachindex(congruent_particles)
            cp = particles(cp_i)
            if p.shape == cp.shape && p.medium == cp.medium
                t_matrices[p_i] = t_matrices[cp_i]
                found = true
                break
            end
        end

        # Congruent particle was not found, we must calculate this t-matrix
        if !found
            t_matrices[p_i] = t_matrix(p.shape, p.medium, medium, ω, Nh)
        end

    end

    return t_matrices
end
