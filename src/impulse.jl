"""
See also: [`ContinuousImpulse`](@ref), [`frequency_to_time`](@ref), [`DiscreteGaussianImpulse`](@ref)

    DiscreteImpulse{T<:AbstractFloat}

A struct used to represent a numerical impulse. Only the fields: `in_freq` which is the frequency response vector, and the frequency vector `ω` are required to use this struct to use in the function `frequency_to_time`.
"""
struct DiscreteImpulse{T<:AbstractFloat}
    t::Vector{T}
    in_time::Vector{Complex{T}}
    ω::Vector{T}
    in_freq::Vector{Complex{T}}
    function DiscreteImpulse{T}(t_vec::AbstractArray{T}, in_time::AbstractVector{Complex{T}}, ω_vec::AbstractArray{T}, in_freq::AbstractVector{Complex{T}}, do_self_test=true) where {T}
        impulse = new{T}(t_vec,collect(in_time),ω_vec,collect(in_freq))
        if do_self_test self_test(impulse) end
        return impulse
    end
end

"""
Check that the discrete impulse vectors are the right sizes
"""
function self_test(impulse::DiscreteImpulse{T}) where {T}

    right_size = length(impulse.t) == length(impulse.in_time)
    right_size = right_size && length(impulse.ω) == length(impulse.in_freq)

    return right_size
end

function DiscreteImpulse(t_vec::AbstractArray{T}, in_time::Union{AbstractArray{Complex{T}},AbstractArray{T}},
        ω_vec::AbstractArray{T} = t_to_ω(t_vec),
        in_freq::Union{AbstractArray{Complex{T}},AbstractArray{T}} = Complex{T}[];
    kws... ) where {T}

    if isempty(in_freq)
        in_freq = collect(time_to_frequency(Vector{T}(in_time), Vector{T}(t_vec), ω_vec; kws...))[:]
    end

    return DiscreteImpulse{T}(
        Vector{T}(t_vec), Vector{Complex{T}}(in_time), Vector{T}(ω_vec), Vector{Complex{T}}(in_freq)
    )
end

"""
See also: [`DiscreteImpulse`](@ref), [`frequency_to_time`](@ref)

    ContinuousImpulse{T<:AbstractFloat}

A struct used to represent an analytic impulse function. Has two fields: `in_time` a function of time `t`, and `in_freq` a function of the angular frequency `ω`. `in_freq` should be the Fourier transform of `in_time`, though this is not enforced.

We use the Fourier transform convention:
F(ω) =  ∫ f(t)*exp(im*ω*t) dt,
f(t) = (2π)^(-1) * ∫ F(ω)*exp(-im*ω*t) dω.

An impluse f(t) is convoluted in time with the field u(t), however we avoid the convolution by working with the fourier transform F(ω) of the impulse f(t), which results in

frequency to time: (2π)^(-1) * ∫ F(ω)*U(ω)*exp(-im*ω*t) dω
"""
struct ContinuousImpulse{T<:AbstractFloat}
    in_time::Function
    in_freq::Function
    # Enforce that the Types are the same
    function ContinuousImpulse{T}(in_time::Function, in_freq::Function, do_self_test=true) where {T}
        impulse = new{T}(in_time,in_freq)
        if do_self_test self_test(impulse) end
        return impulse
    end
end

"""
    continuous_to_discrete_impulse(impulse::ContinuousImpulse, t_vec, ω_vec = t_to_ω(t_vec); t_shift = 0.0, ω_shift = 0.0) 

Returns a [`DiscreteImpulse`](@ref) by sampling `impulse`. The signal can be shifted in time and frequency by choosing `t_shit` and `ω_shift`.
"""
function continuous_to_discrete_impulse(impulse::ContinuousImpulse{T}, t_vec::AbstractArray{T}, ω_vec::AbstractArray{T} = t_to_ω(t_vec); t_shift::T = zero(T), ω_shift::T = zero(T)) where {T}

    # shift in frequency then in time
    in_time = impulse.in_time.(t_vec .- t_shift) .* cos.(ω_shift .* (t_vec .- t_shift))

    # shift in frequency then in time
    in_freq = (impulse.in_freq.(ω_vec .- ω_shift) .+ impulse.in_freq.(ω_vec .+ ω_shift)) .* exp.(im .* ω_vec .* t_shift) ./ T(2)

    return DiscreteImpulse(t_vec, in_time, ω_vec, in_freq)
end

"""
Check that the continuous impulse functions return the correct types
"""
function self_test(impulse::ContinuousImpulse{T}) where {T}

    t = one(T)
    impulse.in_time(t)::Union{T,Complex{T}}

    ω = one(T)
    impulse.in_freq(ω)::Union{T,Complex{T}}

    return true
end

import Base.(+)
function +(i1::ContinuousImpulse{T},i2::ContinuousImpulse{T})::ContinuousImpulse{T} where {T}
    in_time(t) = i1.in_time(t) + i2.in_time(t)
    in_freq(ω) = i1.in_freq(ω) + i2.in_freq(ω)

    ContinuousImpulse{T}(in_time,in_freq)
end

import Base.(*)
function *(a,impulse::ContinuousImpulse{T})::ContinuousImpulse{T} where {T}

    if typeof(one(Complex{T})*a) != Complex{T}
        throw(DomainError(a, "Multiplying impulse by $a would cause impulse type to change, please explicitly cast $a to same type as RegularSource ($T)"))
    end

    in_time(t) = a*impulse.in_time(t)
    in_freq(ω) = a*impulse.in_freq(ω)

    ContinuousImpulse{T}(in_time,in_freq)
end

*(i::ContinuousImpulse,a) = *(a,i)


"""
    TimeDiracImpulse(t0::T)

Dirac Delta function of unit area in the time domain centred at t=t0.

Warning: in the time domain this is a singuarity and so may lead to unexpected behaviour.
"""
function TimeDiracImpulse(t0::T) where {T<:AbstractFloat}
    in_time(t::T) = (t==t0) ? T(Inf) : zero(T)
    in_freq(ω::T) = exp(im*ω*t0)
    ContinuousImpulse{T}(in_time, in_freq, false)
end

function DiscreteTimeDiracImpulse(t0::T, t_vec::AbstractArray{T}, ω_vec::AbstractArray{T} = t_to_ω(t_vec)) where {T}

    return continuous_to_discrete_impulse(TimeDiracImpulse(t0), t_vec, ω_vec)
end

"""
Dirac Delta function of unit area in the frequency domain centred at ω=ω0.

Warning: in frequency space this is a singuarity and so may lead to unexpected behaviour.
"""
function FreqDiracImpulse(ω0::T, dω::T = one(T)) where {T<:AbstractFloat}
    in_time(t::T) = exp(-im*ω0*t) / (T(2)*π)
    in_freq(ω::T) = (ω==ω0) ? T(Inf) : zero(ω)
    ContinuousImpulse{T}(in_time, in_freq, false)
end

"""
See also: [`ContinuousImpulse`](@ref), [`TimeDiracImpulse`](@ref)

    GaussianImpulse(maxω[; σ = 3.0/maxω^2])

Returns a gaussian impulse function, which in the frequency domain is `exp(-σ*ω^2)*(2sqrt(σ*pi))`.
"""
function GaussianImpulse(maxω::T; σ::T = T(3.0)/maxω^2, t_shift = zero(T), ω_shift = zero(T)) where T<:AbstractFloat
    # shift in frequency then in time
    in_time(t::T) = exp(-(t - t_shift)^2 / (4σ)) * cos(ω_shift * (t - t_shift))

    # shift in frequency then in time
    gauss_freq(ω::T) = exp(-σ*ω^2) * (2sqrt(σ*pi))
    in_freq(ω::T) = (gauss_freq(ω - ω_shift) + gauss_freq(ω + ω_shift)) * exp.(im * ω * t_shift) / 2

    ContinuousImpulse{T}(in_time, in_freq)
end

"""
See also: [`ContinuousImpulse`](@ref), [`TimeDiracImpulse`](@ref)

    DiscreteGaussianImpulse(t_vec[, ω_vec])

Returns a discretised gaussian impulse.
"""
function DiscreteGaussianImpulse(t_vec::AbstractArray{T}, ω_vec::AbstractArray{T} = t_to_ω(t_vec);
        kws...) where {T}

    return continuous_to_discrete_impulse(GaussianImpulse(one(T); kws...), t_vec, ω_vec)
end
