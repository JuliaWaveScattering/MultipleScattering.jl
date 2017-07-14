"""
Statistical moments for a batch of models at each value of k or t. Models must
all have the same label (particle radius and volume fraction)
"""
type Moments{T}
    label::Vector{T}
    x_arr::Vector{T}
    moments::Vector{Vector{T}} # First Vector is mean, second is stddev, etc
    num_realisations::Int
end

function Moments{T}(models::Vector{TimeModel{T}}; num_moments=4, response_apply=real)
    # Number of time points sampled
    K = length(models[1].time_arr)
    # Number of realisations
    R = length(models)
    moments = [Vector{T}(K) for i=1:num_moments]
    
    # For each wavenumber, calculate each moment up to num_moments
    for i = 1:K
        responses = [response_apply(models[j].response[i]) for i=1:R]
        μ = mean(responses)
        moment[1][i] = μ
        for m=2:num_moments
            moment = sum((responses .- μ).^m)/(K-1)
            moments[m][i] = sign(moment) * abs(moment)^(1.0/m)
        end
    end
    frequency_models = [models[i].frequency_model for i=1:R]
    label = frequency_models_to_label(frequency_models)
    x_arr = models[1].time_arr
    Moments(label,x_arr,moments,R)
end


function Moments{T}(models::Vector{FrequencyModel{T}}; num_moments=4, response_apply=real)
    # Number of wavenumbers modelled
    K = length(models[1].k_arr)
    # Number of realisations
    R = length(models)
    moments = [Vector{T}(K) for i=1:num_moments]
    
    # For each wavenumber, calculate each moment up to num_moments
    for i = 1:K
        responses = [response_apply(models[j].response[i]) for j=1:R]
        μ = mean(responses)
        moments[1][i] = μ
        for m=2:num_moments
            moment = sum((responses .- μ).^m)/(K-1)
            moments[m][i] = sign(moment) * abs(moment)^(1.0/m)
        end
    end

    label = frequency_models_to_label(models)
    x_arr = models[1].k_arr
    Moments(label,x_arr,moments,R)
end

function frequency_models_to_label{T}(models::Vector{FrequencyModel{T}})
    std_radii = std_radius.(models)
    if mean(std_radii) > 10*eps(T) 
        error("Particles are different sizes: all particles in all models must be the same size to assign the moments a single label.") 
    end
    # We now know they are all the same size, we can just take the first one
    a = models[1].particles[1].r

    volfracs = calculate_volfrac.(models)
    if std(volfracs) > 10*eps(T) 
        error("Models have differen volfracs: all models must have the same volfrac to assign the moments a single label.") 
    end

    return [a,volfracs[1]]
end
