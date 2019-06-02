"""
    statistical_moments(results, n; applytofield=real)::Vector{Matrix}

Calculate moments up to `n` of results at each position and wavenumber/time, after applying `applytofield`.
"""
function statistical_moments(results::AbstractVector{SimRes}, num_moments::Int; applytofield=real) where {T,SimRes<:SimulationResult{T}}

    # Number of positions and wavenumbers/time points sampled
    X, K = size(results[1])

    # Number of realisations
    R = length(results)

    moments = [Matrix{T}(undef,X,K) for i=1:num_moments]

    # For each wavenumber or timestep, calculate each moment up to num_moments
    for i = 1:X
        for j = 1:K
            responses = [applytofield(field(results[r],i,j)) for r=1:R]
            μ = mean(responses)
            moments[1][i,j] = μ
            for m=2:num_moments
                moment = sum((responses .- μ).^m)/(R-1)
                moments[m][i,j] = sign(moment) * abs(moment)^(1.0/m)
            end
        end
    end
    return moments
end
