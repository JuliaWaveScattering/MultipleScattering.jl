"""
See also: [`ContinuousImpulse`](@ref), [`frequency_to_time`](@ref), [`DiscreteGaussianImpulse`](@ref)

    DiscreteImpulse{T<:AbstractFloat}

A struct used to represent a numerical impulse. Only the fields: `in_freq` which is the frequency response vector, and the frequency vector `ω` are required to use this struct to use in the function `frequency_to_time`.
"""
struct DiscreteImpulse{T<:AbstractFloat}
    t::Vector{T}
    in_time::Vector{T}
    ω::Vector{T}
    in_freq::Vector{Complex{T}}
    function DiscreteImpulse{T}(t_vec::AbstractArray{T}, in_time::AbstractVector{T}, ω_vec::AbstractArray{T}, in_freq::AbstractVector{Complex{T}}, do_self_test=true) where {T}
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

function DiscreteImpulse(t_vec::AbstractArray{T}, in_time::AbstractArray{T},
        ω_vec::AbstractArray{T} = t_to_ω(t_vec),
        in_freq::Union{AbstractArray{Complex{T}},Vector{T}} = Complex{T}[];
    kws... ) where {T}

    if isempty(in_freq)
        in_freq = collect(time_to_frequency(Vector{T}(in_time), Vector{T}(t_vec), ω_vec; kws...))[:]
    end

    return DiscreteImpulse{T}(
        Vector{T}(t_vec), Vector{T}(in_time), Vector{T}(ω_vec), Vector{Complex{T}}(in_freq)
    )
end

"""
See also: [`DiscreteImpulse`](@ref), [`frequency_to_time`](@ref)

    ContinuousImpulse{T<:AbstractFloat}

A struct used to represnt analytic impulse functions. Has two fields: `in_time` a function of time `t`, and `in_freq` a function of the angular frequency `ω`. `in_freq` should be the Fourier transform of `in_time`, though this is not enforced.

We use the Fourier transform convention:
F(ω) =  ∫ f(t)*exp(im*ω*t) dt,
f(t) = (2π)^(-1) * ∫ F(ω)*exp(-im*ω*t) dω.

An impluse f(t) is convoluted in time with the field u(t), however we avoid the convlution by working with the fourier transform F(ω) of the impulse f(t), which results in

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

function continuous_to_discrete_impulse(impulse::ContinuousImpulse{T}, t_vec::AbstractArray{T}, ω_vec::AbstractArray{T} = t_to_ω(t_vec)) where {T}

    return DiscreteImpulse(t_vec, impulse.in_time.(t_vec), ω_vec, impulse.in_freq.(ω_vec))
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
        throw(DomainError(a, "Multiplying impulse by $a would cause impulse type to change, please explicitly cast $a to same type as Source ($T)"))
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

    GaussianImpulse(maxω[, a = 3.0/maxω^2])

Returns a gaussian impulse function, which in the frequency domain is exp(-a*ω^2)*(2sqrt(a*pi)).
"""
function GaussianImpulse(maxω::T, a::T = T(3.0)/maxω^2; time_shift = zero(T)) where {T<:AbstractFloat}
    in_time(t::T) = exp(-(t-time_shift)^2/(4a))
    in_freq(ω::T) = exp(-a*ω^2)*exp(im*ω*time_shift)*(2sqrt(a*pi))
    ContinuousImpulse{T}(in_time, in_freq)
end

"""
See also: [`ContinuousImpulse`](@ref), [`TimeDiracImpulse`](@ref)

    DiscreteGaussianImpulse(t_vec[, ω_vec])

Returns a discretised gaussian impulse.
"""
function DiscreteGaussianImpulse(t_vec::AbstractArray{T}, ω_vec::AbstractArray{T} = t_to_ω(t_vec);
        a::T = T(2.48)/maximum(ω_vec)^2) where {T}

    return continuous_to_discrete_impulse(GaussianImpulse(one(T), a), t_vec, ω_vec)
end
