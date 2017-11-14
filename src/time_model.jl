
type TimeModel{T}
    frequency_model::FrequencyModel{T}
    response::Matrix{T}
    time_arr::Vector{T}
    impulse::Function
end

"Convert a frequency model into a time model using the inverse fourier transform"
function TimeModel{T}(
        freq_model::FrequencyModel{T};
        timesteps=length(freq_model.k_arr),
        time_arr=linspace(zero(T),one(T)*π*2/firstnonzero(freq_model.k_arr),timesteps+1)[1:timesteps],
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
Calculates the time response from the frequency response by approximating an inverse Fourier transform.
 The time signal is assumed to be real and only positive frequenices k_arr given.
 The result is convoluted with the user specified impulse, which is a function of the frequency.
"""
function frequency_to_time{T}(freq_response::Matrix{Complex{T}},k_arr::AbstractArray{T},
  time_arr::AbstractArray{T}, impulse::Function;
  addzerofrequency=true)
    # we assume the convention: f(t) = 1/(2π) ∫f(w)exp(-im*w*t)dw

    if addzerofrequency && !(k_arr[1] ≈ 0) # adds the response for k=0 if not present
      f0 = (k_arr[2]*freq_response[1] - k_arr[1]*freq_response[2])/(k_arr[2]-k_arr[1])
      freq_response = [f0; freq_response]
      k_arr = [0; k_arr]
    end
    timesteps = length(time_arr)
    positions = size(freq_response,2)
    freqsteps = size(freq_response,1)

    impulse_vec = map(impulse,k_arr)

    u = Matrix{Complex{T}}(timesteps,positions)
    for i=1:timesteps
        for j=1:positions
            u[i,j] = zero(T)
            # using trapezoidal rule
            for ki=1:(freqsteps-1)
                dk = k_arr[ki+1] - k_arr[ki]
                k = k_arr[ki]/2 + k_arr[ki+1]/2
                uhat = (freq_response[ki]+freq_response[ki+1])/2
                t = time_arr[i]
                u[i,j] += impulse(k)*uhat*exp(-im*k*t)*dk
            end
        end
    end

    return real(u)/pi # constant due to our convention and because we are using only positive frequencies.
end
