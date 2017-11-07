
type TimeModel{T}
    frequency_model::FrequencyModel{T}
    response::Matrix{Complex{T}}
    time_arr::Vector{T}
    impulse::Function
end

"Convert a frequency model into a time model using the inverse fourier transform"
function TimeModel{T}(
        freq_model::FrequencyModel{T};
        timesteps=length(freq_model.k_arr)+1,
        time_arr=linspace(zero(T),one(T)*π*2/firstnonzero(freq_model.k_arr),timesteps),
        impulse=delta_fnc
    )
    response = frequency_to_time(freq_model.response,freq_model.k_arr,time_arr,impulse)
    TimeModel{T}(freq_model,response,time_arr,impulse)
end


"""
Returns the first element of array which isn't zero (assumes elements are 
increasing and distinct)
"""
function firstnonzero{T}(arr::AbstractArray{T})
    if arr[1] != 0
        return arr[1]
    else
        return arr[2]
    end
end


"""
Function which is one everywhere in frequency domain. Represents a delta 
function of unit area in the time domain, centred at t=zero.
"""
delta_fnc{T}(k::T) = one(T)


"""
Perform an inverse Fourier transform with the frequency response.
Calculates the response of the system to the frequency impulse, used to convert
a frequency model into a time model.
"""
function frequency_to_time{T}(freq_response::Matrix{Complex{T}},k_arr::AbstractArray{T},time_arr::AbstractArray{T},impulse::Function)

    timesteps = length(time_arr)
    positions = size(freq_response,2)
    freqsteps = size(freq_response,1)

    impulse_vec = map(impulse,k_arr)

    # Size of each slice of k which we integrate over, assumes k_arr is uniform
    dk = k_arr[2] - k_arr[1]

    u = Matrix{Complex{T}}(timesteps,positions)
    for i=1:timesteps
        for j=1:positions
            u[i,j] = zero(T)
            for ki=1:freqsteps
                k = k_arr[ki]
                uhat = freq_response[ki]
                t = time_arr[i]
                u[i,j] += impulse(k)*uhat*exp(im*k*t)*dk
            end
        end
    end

    return u
end


