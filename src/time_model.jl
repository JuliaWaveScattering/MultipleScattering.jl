
type TimeModel{T}
    simulation::FrequencyModel{T}
    response::Matrix{Complex{T}}
    time::Vector{T}
end

"Convert a frequency model into a time model using the inverse fourier transform"
function TimeModel{T}(freq_model::FrequencyModel{T})
    warn("FrequencyModel to TimeModel converter not implemented yet")
end