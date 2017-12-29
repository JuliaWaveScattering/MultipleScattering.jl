type TimeSimulation{T}
    frequency_simulation::FrequencySimulation{T}
    response::Matrix{T}
    time_arr::Vector{T}
    impulse_arr::Vector{Complex{T}}
    method::Symbol
end

"""
Convert a frequency simulation into a time simulation using the inverse fourier transform.
Assumes only positive frequencies and a real time signal
"""
function TimeSimulation{T}(
        freq_simulation::FrequencySimulation{T};
        time_arr = ω_to_t(freq_simulation.k_arr),
        impulse = get_gaussian_freq_impulse(maximum(freq_simulation.k_arr)),
        method =:dft
    )
    response = frequency_to_time(freq_simulation.response,freq_simulation.k_arr,time_arr,impulse; method = method)
    TimeSimulation{T}(freq_simulation,response,time_arr,impulse.(freq_simulation.k_arr),method)
end

"Take simulation parameters, run simulation and populate the response array."
function generate_responses!{T}(TimeSimulation::TimeSimulation{T},t_arr::Vector{T}=TimeSimulation.time_arr)
    TimeSimulation.time_arr = t_arr
    TimeSimulation.response = frequency_to_time(
        TimeSimulation.frequency_simulation.response, TimeSimulation.frequency_simulation.k_arr,
        TimeSimulation.time_arr, TimeSimulation.impulse; method = TimeSimulation.method
    )
end

"""
Calcuate the volume fraction of this simulation from the volume of the particles
and the bounding shape.
"""
function calculate_volfrac(time_simulation::TimeSimulation)
    volume(time_simulation.frequency_simulation.particles)/volume(time_simulation.frequency_simulation.shape)
end


"""
returns an array of time from the frequency array ω_arr.
Uses the same convention for sampling the time as the discrete Fourier transfrom.
Assumes ω_arr is ordered and non-negative.
"""
function ω_to_t{T}(ω_arr::AbstractArray{T})
  N = length(ω_arr)
  if (ω_arr[1] == 0.)
    N -= 1
  elseif minimum(ω_arr)<0
    error("expected only non-negative values for the frequencies")
  end
  dω = ω_arr[2]-ω_arr[1]
  time_arr = linspace(0,2π/dω,2N+2)[1:(2N+1)]
  return time_arr
end

"The inverse of ω_to_t if ω_arr[1] == 0"
function t_to_ω{T}(t_arr::AbstractArray{T})
  N = Int((length(t_arr)-1)/2)
  maxt = t_arr[2]-t_arr[1] + t_arr[end]
  maxω = N*2π/maxt
  ω_arr = linspace(0,maxω,N+1)
  return ω_arr
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
delta_freq_impulse{T}(k::T) = one(T)

"""
Returns a gaussian impulse function in the frequency domain. In the time domain
this impulse is exp(-t^2/(4a))
"""
get_gaussian_freq_impulse{T}(maxk::T, a::T = T(2.48)/maxk^2) = ω -> exp(-a*ω^2)*(2sqrt(a*pi))

"""
Returns a gaussian impulse function in the time domain.
"""
get_gaussian_time_impulse{T}(maxk::T, a::T = T(2.48)/maxk^2) = t -> exp(-t^2/(4a))


"""
Calculates the time response from the frequency response by approximating an
inverse Fourier transform. The time signal is assumed to be real and the
frequenices k_arr are assumed to be positive (can include zero) and sorted. The
result is convoluted ωith the user specified impulse, which is a function of the
frequency.
"""
function frequency_to_time{T}(freq_response::Matrix{Complex{T}},ω_arr::AbstractArray{T},
  time_arr::AbstractArray{T}=ω_to_t(ω_arr), impulse::Function = delta_freq_impulse;
  addzerofrequency=true, method=:dft)
    # we assume the convention: f(t) = 1/(2π) ∫f(ω)exp(-im*ω*t)dω

    positions = size(freq_response,2)

    # if k=0 is not present, then add the response for k=0 by using linear interpolation.
    if addzerofrequency && minimum(ω_arr) > 0.0
      zeroresponse = [
        (ω_arr[2]*freq_response[1,j] - ω_arr[1]*freq_response[2,j])/(ω_arr[2]-ω_arr[1])
      for j in 1:positions]
      freq_response = [transpose(zeroresponse); freq_response]
      ω_arr = [0; ω_arr]
    end

    # determine how to approximate  ∫f(ω)exp(-im*ω*t)dω
    freqsteps = size(freq_response,1)
    if method == :trapezoidal
      fourier_integral = (t,j) ->
        sum(1:(freqsteps-1)) do ωi
            ω = ω_arr[ωi]
            uhat = freq_response[ωi,j]
            f1 = impulse(ω)*uhat*exp(-im*ω*t)
            ω = ω_arr[ωi+1]
            uhat = freq_response[ωi+1,j]
            f2 = impulse(ω)*uhat*exp(-im*ω*t)
            dω = ω_arr[ωi+1] - ω_arr[ωi]
            (f1+f2)*dω/2
        end
    else
      # integral to match exactly standard dft
      fourier_integral = (t,j) ->
        impulse(0)*freq_response[1,j]*(ω_arr[2]-ω_arr[1])/2 +  # because ω=0 should not be added tωice later
        sum(2:freqsteps) do ωi
          dω = ω_arr[ωi]-ω_arr[ωi-1]
          ω = ω_arr[ωi]
          uhat = freq_response[ωi,j]
          impulse(ω)*uhat*exp(-im*ω*t)*dω
        end
    end
    # impulse(ω_arr[freqsteps])*freq_response[freqsteps,j]*(ω_arr[freqsteps]-ω_arr[freqsteps-1]) +

    u = [ fourier_integral(t,j) for t in time_arr, j=1:positions]

    return real(u)/pi # constant due to our Fourier convention and because only positive frequencies are used.
end

"""
The inverse of the function frequency_to_time (only an exact inverse when using
:dft integration)
"""
function time_to_frequency{T}(time_response::Matrix{T},ω_arr::AbstractArray{T},
  time_arr::AbstractArray{T}=ω_to_t(ω_arr); impulse::Function = delta_freq_impulse, method=:dft)
    # we assume the convention: f(ω) =  ∫f(t)exp(im*ω*t)dt

    # determine how to approximate  ∫f(t)exp(im*ω*t)dt
    timesteps = size(time_response,1)
    if method == :trapezoidal
      fourier_integral = (ω,j) ->
        sum(1:(timesteps-1)) do ti
            dt = time_arr[ti+1] - time_arr[ti]
            t = time_arr[ti]
            u = time_response[ti,j]
            f1 = impulse(t)*u*exp(im*ω*t)
            t = time_arr[ti+1]
            u = time_response[ti+1,j]
            f2 = impulse(t)*u*exp(im*ω*t)
            (f1+f2)*dt/2
        end
    else
      # integral to match exactly standard dft
      fourier_integral = (ω,j) ->
        impulse(0)*time_response[1,j]*(time_arr[2]-time_arr[1]) +
        sum(2:timesteps) do ti
          dt = time_arr[ti]-time_arr[ti-1]
          t  = time_arr[ti]
          u  = time_response[ti,j]
          impulse(t)*u*exp(im*ω*t)*dt
        end
    end

    positions = size(time_response,2)
    uhat = [fourier_integral(ω,j) for ω in ω_arr, j=1:positions]

    return uhat
end
