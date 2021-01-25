
"Create the matrix S which will be inverted to find the scattering coefficients."
function scattering_matrix(medium::PhysicalMedium, particles::AbstractParticles, t_matrices::Vector, ω::T, order::Integer)::Matrix{Complex{T}} where T
    # Generate response for one specific k
    # Number of particles
    P = length(particles)

    # Length of scattering basisfor each particle
    N = basisorder_to_basislength(typeof(medium),order)

    # No particles means no scattering
    if P == 0
        return Matrix{Complex{T}}(undef,0,0)
    end

    # Faire: this could potentially return an MMatrix
    function S_block(j,l)
        if j == l
            return zeros(Complex{T}, N, N)
        else
            x_lj = origin(particles[j]) .- origin(particles[l])
            U = outgoing_translation_matrix(medium, order, ω, x_lj)
            return - transpose(U) * t_matrices[l]
        end
    end

    # Evaluates all the blocks in the big matrix
    S_blocks = [S_block(j,l) for j in 1:P, l in 1:P]

    # Reshape S_blocks into big matrix
    S = zeros(Complex{T}, P*N, P*N)
    for i in 1:P
        for j in 1:P
            S[((i-1)*N+1):(i*N), ((j-1)*N+1):(j*N)] .= S_blocks[i,j]
        end
    end

    return S
end
