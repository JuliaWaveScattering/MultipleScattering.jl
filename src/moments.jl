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

function StatisticalMoments{T}(models::Vector{TimeSimulation{T}}, num_moments=4; response_apply=real)
    StatisticalMoments(
        models_to_label(models),
        models[1].time_arr,
        calculate_moments(models,num_moments),
        length(models)
    )
end

function StatisticalMoments{T}(models::Vector{FrequencySimulation{T}}, num_moments=4; response_apply=real)
    StatisticalMoments(
        models_to_label(models),
        models[1].k_arr,
        calculate_moments(models,num_moments;response_apply=response_apply),
        length(models)
    )
end

function calculate_moments(models::Vector,num_moments::Int;response_apply=real)
    # Number of wavenumbers or time points sampled
    K = length(models[1].response)
    # Number of realisations
    R = length(models)
    moments = [Vector{typeof(models[1].response[1,1].re)}(K) for i=1:num_moments]

    # For each wavenumber or timestep, calculate each moment up to num_moments
    for i = 1:K
        responses = [response_apply(models[j].response[i,1]) for j=1:R]
        μ = mean(responses)
        moments[1][i] = μ
        for m=2:num_moments
            moment = sum((responses .- μ).^m)/(R-1)
            moments[m][i] = sign(moment) * abs(moment)^(1.0/m)
        end
    end
    return moments
end

function models_to_label{T}(models::Vector{TimeSimulation{T}})
    frequency_models = [models[i].frequency_model for i=1:R]
    return models_to_label(frequency_models)
end

function models_to_label{T}(models::Vector{FrequencySimulation{T}})
    std_radii = std_radius.(models)
    if mean(std_radii) > 10*eps(T)
        error("Particles are different sizes: all particles in all models must be the same size to assign the moments a single label.")
    end
    # We now know they are all the same size, we can just take the first one
    a = models[1].particles[1].r

    volfracs = calculate_volfrac.(models)
    if std(volfracs) > 10*eps(T)
        error("Models have different volfracs: all models must have the same volfrac to assign the moments a single label.")
    end

    return [a,volfracs[1]]
end
