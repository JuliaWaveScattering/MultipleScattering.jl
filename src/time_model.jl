
type TimeModel{T}
    frequency_model::FrequencyModel{T}
    response::Matrix{T}
    time_arr::Vector{T}
    impulse::Function
    method::Symbol
end

"Convert a frequency model into a time model using the inverse fourier transform. Assumes only positive frequencies and a real time signal"
function TimeModel{T}(
        freq_model::FrequencyModel{T};
        time_arr = wTot(freq_model.k_arr),
        impulse = gaussian_impulses(maximum(freq_model.k_arr)),
        method =:dft
    )
    response = frequency_to_time(freq_model.response,freq_model.k_arr,time_arr,impulse; method = method)
    TimeModel{T}(freq_model,response,time_arr,impulse,method)
end

"Take model parameters, run model and populate the response array."
function generate_responses!{T}(timemodel::TimeModel{T},t_arr::Vector{T}=timemodel.time_arr)
    timemodel.time_arr = t_arr
    timemodel.response = frequency_to_time(
        timemodel.frequency_model.response, timemodel.frequency_model.k_arr,
        timemodel.time_arr, timemodel.impulse; method = timemodel.method
    )
end

"""
returns an array of time from the frequency array w_arr.
Uses the same convention for sampling the time as the discrete Fourier transfrom.
Assumes w_arr is ordered and non-negative.
"""
function wTot{T}(w_arr::AbstractArray{T})
  N = length(w_arr)
  if (w_arr[1] == 0.)
    N -= 1
  elseif minimum(w_arr)<0
    error("expected only non-negative values for the frequencies")
  end
  dw = w_arr[2]-w_arr[1]
  time_arr = linspace(0,2π/dw,2N+2)[1:(2N+1)]
  return time_arr
end

"The inverse of wTot if w_arr[1] == 0"
function tTow{T}(t_arr::AbstractArray{T})
  N = Int((length(t_arr)-1)/2)
  maxt = t_arr[2]-t_arr[1] + t_arr[end]
  maxw = N*2π/maxt
  w_arr = linspace(0,maxw,N+1)
  return w_arr
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
Returns a gaussian impulse function in the frequency domain. In the time domain this impulse is exp(-t^2/(4a))
"""
gaussian_impulses{T}(maxk::T, a::T = T(2.48)/maxk^2) = w -> exp(-a*w^2)*(2sqrt(a*pi))


"""
Calculates the time response from the frequency response by approximating an inverse Fourier transform.
 The time signal is assumed to be real and the frequenices k_arr are assumed to be positive (can include zero) and sorted.
 The result is convoluted with the user specified impulse, which is a function of the frequency.
"""
function frequency_to_time{T}(freq_response::Matrix{Complex{T}},w_arr::AbstractArray{T},
  time_arr::AbstractArray{T}=wTot(w_arr), impulse::Function = delta_fnc;
  addzerofrequency=true, method=:dft)
    # we assume the convention: f(t) = 1/(2π) ∫f(w)exp(-im*w*t)dw

    positions = size(freq_response,2)

    # if k=0 is not present, then add the response for k=0 by using linear interpolation.
    if addzerofrequency && minimum(w_arr) > 0.0
      zeroresponse = [
        (w_arr[2]*freq_response[1,j] - w_arr[1]*freq_response[2,j])/(w_arr[2]-w_arr[1])
      for j in 1:positions]
      freq_response = [transpose(zeroresponse); freq_response]
      w_arr = [0; w_arr]
    end

    # determine how to approximate  ∫f(w)exp(-im*w*t)dw
    freqsteps = size(freq_response,1)
    if method == :trapezoidal
      fourier_integral = (t,j) ->
        sum(1:(freqsteps-1)) do wi
            w = w_arr[wi]
            uhat = freq_response[wi,j]
            f1 = impulse(w)*uhat*exp(-im*w*t)
            w = w_arr[wi+1]
            uhat = freq_response[wi+1,j]
            f2 = impulse(w)*uhat*exp(-im*w*t)
            dw = w_arr[wi+1] - w_arr[wi]
            (f1+f2)*dw/2
        end
    else
      # integral to match exactly standard dft
      fourier_integral = (t,j) ->
        impulse(0)*freq_response[1,j]*(w_arr[2]-w_arr[1])/2 +  # because w=0 should not be added twice later
        sum(2:freqsteps) do wi
          dw = w_arr[wi]-w_arr[wi-1]
          w = w_arr[wi]
          uhat = freq_response[wi,j]
          impulse(w)*uhat*exp(-im*w*t)*dw
        end
    end
    # impulse(w_arr[freqsteps])*freq_response[freqsteps,j]*(w_arr[freqsteps]-w_arr[freqsteps-1]) +

    u = [ fourier_integral(t,j) for t in time_arr, j=1:positions]

    return real(u)/pi # constant due to our Fourier convention and because only positive frequencies are used.
end

"""
The inverse of the function frequency_to_time (only an exact inverse when using :riemann integration)
"""
function time_to_frequency{T}(time_response::Matrix{T},w_arr::AbstractArray{T},
  time_arr::AbstractArray{T}=wTot(w_arr); impulse::Function = delta_fnc, method=:dft)
    # we assume the convention: f(w) =  ∫f(t)exp(im*w*t)dt

    # determine how to approximate  ∫f(t)exp(im*w*t)dt
    timesteps = size(time_response,1)
    if method == :trapezoidal
      fourier_integral = (w,j) ->
        sum(1:(timesteps-1)) do ti
            dt = time_arr[ti+1] - time_arr[ti]
            t = time_arr[ti]
            u = time_response[ti,j]
            f1 = impulse(t)*u*exp(im*w*t)
            t = time_arr[ti+1]
            u = time_response[ti+1,j]
            f2 = impulse(t)*u*exp(im*w*t)
            (f1+f2)*dt/2
        end
    else
      # integral to match exactly standard dft
      fourier_integral = (w,j) ->
        impulse(0)*time_response[1,j]*(time_arr[2]-time_arr[1]) +
        sum(2:timesteps) do ti
          dt = time_arr[ti]-time_arr[ti-1]
          t  = time_arr[ti]
          u  = time_response[ti,j]
          impulse(t)*u*exp(im*w*t)*dt
        end
    end

    positions = size(time_response,2)
    uhat = [fourier_integral(w,j) for w in w_arr, j=1:positions]

    return uhat # constant due to our convention and because we are using only positive frequencies.
end
