
# mutable struct TimeSimulation{Dim,P<:PhysicalProperties,T<:AbstractFloat} <: Simulation{Dim,T}
#     medium::P
#     particles::Particles{T,Dim}
#     impulse::Function
#     source::Source{P,T}
# end

# Main run function, all other run functions use this
function TimeSimulationResult(simres::FrequencySimulationResult{T,Dim,FieldDim};
        t_vec = ω_to_t(simres.ω),
        impulse = get_gaussian_freq_impulse(maximum(simres.ω)),
        impulse_vec = impulse.(simres.ω),
        method =:dft
    ) where {Dim,FieldDim,T}

    time_field = frequency_to_time(field(simres), simres.ω, t_vec;
        impulse = impulse, impulse_vec = impulse_vec, method = method)

    return TimeSimulationResult(time_field, simres.x, t_vec)
end

#
# """
# Convert a frequency simulation into a time simulation using the inverse fourier transform.
# Assumes only positive frequencies and a real time signal
# """
# function TimeSimulation(
#         freq_simulation::FrequencySimulation{T};
#         t_vec = ω_to_t(freq_simulation.k_arr),
#         impulse = get_gaussian_freq_impulse(maximum(freq_simulation.k_arr)),
#         impulse_vec = impulse.(freq_simulation.k_arr),
#         method =:dft
#     ) where T <: AbstractFloat
#     response = frequency_to_time(freq_simulation.response,freq_simulation.k_arr,t_vec;
#         impulse = impulse, impulse_vec = impulse_vec, method = method)
#
#     real(u)/pi
#
#     TimeSimulation{T}(freq_simulation,response,t_vec,impulse_vec,method)
# end
#
# "Take simulation parameters, run simulation and populate the response array."
# function generate_responses!{T}(TimeSimulation::TimeSimulation{T},t_arr::Vector{T}=TimeSimulation.t_vec)
#     TimeSimulation.t_vec = t_arr
#     TimeSimulation.response = frequency_to_time(
#         TimeSimulation.frequency_simulation.response, TimeSimulation.frequency_simulation.k_arr,
#         TimeSimulation.t_vec; impulse_vec = TimeSimulation.impulse_vec, method = TimeSimulation.method
#     )
# end

"""
returns an array of time from the frequency array ω_vec.
Uses the same convention for sampling the time as the discrete Fourier transfrom.
Assumes ω_vec is ordered and non-negative.
"""
function ω_to_t(ω_vec::AbstractArray{T}) where T <: AbstractFloat
  N = length(ω_vec)
  if (ω_vec[1] == zero(T))
    N -= 1
  elseif minimum(ω_vec) < zero(T)
    error("expected only non-negative values for the frequencies")
  end
  dω = ω_vec[2] - ω_vec[1]
  t_vec = linspace(zero(T),2π/dω,2N+2)[1:(2N+1)]
  return t_vec
end

"The inverse of ω_to_t if ω_vec[1] == 0"
function t_to_ω(t_arr::AbstractVector{T}) where T <: AbstractFloat
  N = Int((length(t_arr)-one(T))/2)
  maxt = t_arr[2]-t_arr[1] + t_arr[end]
  maxω = N*2π/maxt
  ω_vec = linspace(zero(T),maxω,N+1)
  return ω_vec
end

"""
Returns the first element of array which isn't zero (assumes elements are
increasing and distinct)
"""
function firstnonzero(arr::AbstractArray{T}) where T <: AbstractFloat
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
delta_freq_impulse(k::T) where T <: AbstractFloat = one(T)

"""
Returns a gaussian impulse function in the frequency domain. In the time domain
this impulse is exp(-t^2/(4a))
"""
get_gaussian_freq_impulse(maxk::T, a::T = T(2.48)/maxk^2) where T <: AbstractFloat = ω -> exp(-a*ω^2)*(2sqrt(a*pi))

"""
Returns a gaussian impulse function in the time domain.
"""
get_gaussian_time_impulse(maxk::T, a::T = T(2.48)/maxk^2) where T <: AbstractFloat = t -> exp(-t^2/(4a))


"""
Calculates the time response from the frequency response by approximating an
inverse Fourier transform. The time signal is assumed to be real and the
frequenices k_arr are assumed to be positive (can include zero) and sorted. The
result is convoluted ωith the user specified impulse, which is a function of the
frequency.
"""
function frequency_to_time(field_mat::Matrix{Complex{T}}, ω_vec::Union{AbstractVector{T},RowVector{T}},
        t_vec::Union{AbstractVector{T},RowVector{T}} = ω_to_t(ω_vec);
        impulse::Function = delta_freq_impulse,
        impulse_vec = impulse.(ω_vec),
        addzerofrequency=true, method=:dft) where T <: AbstractFloat
    # we assume the convention: f(t) = 1/(2π) ∫f(ω)exp(-im*ω*t)dω

    positions = size(field_mat,1)

    # if k=0 is not present, then add the response for k=0 by using linear interpolation.
    if addzerofrequency && minimum(ω_vec) > zero(T)
      zeroresponse = [
        (ω_vec[2]*field_mat[j,1] - ω_vec[1]*field_mat[j,2])/(ω_vec[2]-ω_vec[1])
      for j in 1:positions]

      field_mat = reshape([zeroresponse; field_mat[:]], (positions, size(field_mat,2)+1) )

      if impulse_vec == impulse.(ω_vec)
        impulse_vec = impulse.([0; ω_vec])
      else
        impulse_zero = (ω_vec[2]*impulse_vec[1] - ω_vec[1]*impulse_vec[2])/(ω_vec[2]-ω_vec[1])
        impulse_vec = impulse.([impulse_zero; impulse_vec])
      end
      ω_vec = [0; ω_vec]
    end

    # determine how to approximate  ∫f(ω)exp(-im*ω*t)dω
    ω_steps = size(field_mat,2)
    if method == :trapezoidal
      fourier_integral = (t,j) ->
        sum(1:(ω_steps-1)) do ωi
            ω = ω_vec[ωi]
            uhat = field_mat[j,ωi]
            f1 = impulse_vec[ωi]*uhat*exp(-im*ω*t)
            ω = ω_vec[ωi+1]
            uhat = field_mat[j,ωi+1]
            f2 = impulse_vec[ωi+1]*uhat*exp(-im*ω*t)
            dω = ω_vec[ωi+1] - ω_vec[ωi]
            (f1+f2)*dω/2
        end
    else
      # integral to match exactly standard dft, assuming ω_vec[1] == 0
      fourier_integral = (t,j) ->
        impulse_vec[1]*field_mat[j,1]*(ω_vec[2]-ω_vec[1])/2 +  # because ω=0 should not be added tωice later
        sum(2:ω_steps) do ωi
          dω = ω_vec[ωi]-ω_vec[ωi-1]
          ω = ω_vec[ωi]
          uhat = field_mat[j,ωi]
          impulse_vec[ωi]*uhat*exp(-im*ω*t)*dω
        end
    end
    # impulse(ω_vec[ω_steps])*field_mat[ω_steps,j]*(ω_vec[ω_steps]-ω_vec[ω_steps-1]) +

    u = [ fourier_integral(t,j) for j=1:positions, t in t_vec]

    return real.(u)/pi # constant due to our Fourier convention and because only positive frequencies are used.
end

"""
The inverse of the function frequency_to_time (only an exact inverse when using
:dft integration)
"""
function time_to_frequency(time_response::Matrix{T},ω_vec::AbstractArray{T},
  t_vec::AbstractArray{T} = ω_to_t(ω_vec);
  impulse::Function = delta_freq_impulse,
  impulse_vec = impulse.(t_vec),
  method=:dft) where T <: AbstractFloat
    # we assume the convention: f(ω) =  ∫f(t)exp(im*ω*t)dt

    # determine how to approximate  ∫f(t)exp(im*ω*t)dt
    timesteps = size(time_response,1)
    if method == :trapezoidal
      fourier_integral = (ω,j) ->
        sum(1:(timesteps-1)) do ti
            dt = t_vec[ti+1] - t_vec[ti]
            t = t_vec[ti]
            u = time_response[ti,j]
            f1 = impulse_vec[ti]*u*exp(im*ω*t)
            t = t_vec[ti+1]
            u = time_response[ti+1,j]
            f2 = impulse_vec[ti+1]*u*exp(im*ω*t)
            (f1+f2)*dt/2
        end
    else
      # integral to match exactly standard dft
      fourier_integral = (ω,j) ->
        impulse(0)*time_response[1,j]*(t_vec[2]-t_vec[1]) +
        sum(2:timesteps) do ti
          dt = t_vec[ti]-t_vec[ti-1]
          t  = t_vec[ti]
          u  = time_response[ti,j]
          impulse_vec[ti]*u*exp(im*ω*t)*dt
        end
    end

    positions = size(time_response,2)
    uhat = [fourier_integral(ω,j) for ω in ω_vec, j=1:positions]

    return uhat
end
