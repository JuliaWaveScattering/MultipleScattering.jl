
type Moments{T}
    label::Array{T}
    x_arr::Array{T}
    moments::Array{Array{T}} # Array of moments, first array is the mean, second the variance...
    num_realisations::UInt
end

function Moments{T}(models::Array{TimeModel{T}}; num_moments=4, resp_apply=real)
    resps = [resp_apply(model.response) for model in models][:]
    m = mean(resps)[:]
    L = length(resps) - 1
    moments = [ sum((r[:]-m).^n for r in resps) ./ L for  n = 2:num_moments ]
    moments = [ sign(moments[i]) .* abs(moments[i]).^(1 / (i + 1)) for i in indices(moments, 1)]
    # moments[1] = sqrt(moments[1])
    # moments[3:end] = [ moments[i]./(moments[1].^(i+2.)) for i in indices(moments[3:end],1)]
    Moments(modelTolabel(models[1]), models[1].x_arr, [m, moments...], length(models)) # assumes each model in models has the same label
end

function Moments{T}(models::Array{FrequencyModel{T}}; num_moments=4, resp_apply=real)
    resps = [resp_apply(model.response) for model in models][:]
    m = mean(resps)[:]
    L = length(resps) - 1
    moments = [ sum((r[:]-m).^n for r in resps) ./ L for  n = 2:num_moments ]
    moments = [ sign(moments[i]) .* abs(moments[i]).^(1 / (i + 1)) for i in indices(moments,1)]
    # moments[1] = sqrt(moments[1])
    # moments[3:end] = [ moments[i]./(moments[1].^(i+2.)) for i in indices(moments[3:end],1)]
    Moments(modelTolabel(models[1]), models[1].x_arr, [m, moments...], length(models)) # assumes each model in models has the same label
end

function Moments(label::Array{Float64},moments::Array{Array{Float64,1},1})
    Moments(label, 1.0:1.0:length(moments[1]), moments, 0)
end