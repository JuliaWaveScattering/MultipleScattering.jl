"""
Statistical moments for a batch of simulations at each value of k or t. Simulations must
all have the same label (particle radius and volume fraction)
"""
type StatisticalMoments{T}
    label::Vector{T}
    x_arr::Vector{T}
    moments::Vector{Vector{T}} # First Vector is mean, second is stddev, etc
    num_realisations::Int
end

function StatisticalMoments{T}(simulations::Vector{TimeSimulation{T}}, num_moments=4; response_apply=real)
    StatisticalMoments(
        simulations_to_label(simulations),
        simulations[1].time_arr,
        calculate_moments(simulations,num_moments),
        length(simulations)
    )
end

function StatisticalMoments{T}(simulations::Vector{FrequencySimulation{T}}, num_moments=4; response_apply=real)
    StatisticalMoments(
        simulations_to_label(simulations),
        simulations[1].k_arr,
        calculate_moments(simulations,num_moments;response_apply=response_apply),
        length(simulations)
    )
end

function calculate_moments(simulations::Vector,num_moments::Int;response_apply=real)
    # Number of wavenumbers or time points sampled
    K = length(simulations[1].response)
    # Number of realisations
    R = length(simulations)
    moments = [Vector{typeof(response_apply(simulations[1].response[1,1]))}(K) for i=1:num_moments]

    # For each wavenumber or timestep, calculate each moment up to num_moments
    for i = 1:K
        responses = [response_apply(simulations[j].response[i,1]) for j=1:R]
        μ = mean(responses)
        moments[1][i] = μ
        for m=2:num_moments
            moment = sum((responses .- μ).^m)/(R-1)
            moments[m][i] = sign(moment) * abs(moment)^(1.0/m)
        end
    end
    return moments
end

function simulations_to_label{T}(simulations::Vector{TimeSimulation{T}})
    frequency_simulations = [s.frequency_simulation for s in simulations]
    return simulations_to_label(frequency_simulations)
end

function simulations_to_label{T}(simulations::Vector{FrequencySimulation{T}})
    std_radii = std_radius.(simulations)
    if mean(std_radii) > 100*eps(T)
        error("Particles are different sizes: all particles in all simulations must be the same size to assign the moments a single label.")
    end
    # We now know they are all the same size, we can just take the first one
    a = simulations[1].particles[1].r

    # volfracs are never exactly the same
    volfracs = calculate_volfrac.(simulations)
    if std(volfracs) > T(0.01)
        error("Simulations have different volfracs: all simulations must have the same volfrac to assign the moments a single label.")
    end

    return [mean(volfracs),a]
end
