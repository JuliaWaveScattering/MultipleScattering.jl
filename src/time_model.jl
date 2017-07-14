
type TimeModel{T}
    frequency_model::FrequencyModel{T}
    response::Matrix{Complex{T}}
    time_arr::Vector{T}
end

"Convert a frequency model into a time model using the inverse fourier transform"
function TimeModel{T}(freq_model::FrequencyModel{T})
    warn("FrequencyModel to TimeModel converter not implemented yet")
end