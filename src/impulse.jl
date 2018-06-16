
struct Impulse{T<:AbstractFloat}
    in_time::Function
    in_freq::Function
    # Enforce that the Types are the same
    function Impulse{T}(in_time::Function, in_freq::Function, do_self_test=true) where {T}
        impulse = new{T}(in_time,in_freq)
        if do_self_test self_test(impulse) end
        return impulse
    end
end

"""
Check that the impulse functions return the correct types
"""
function self_test(impulse::Impulse{T}) where {T}

    t = one(T)
    impulse.in_time(t)::Union{T,Complex{T}}

    ω = one(T)
    impulse.in_freq(ω)::Union{T,Complex{T}}

    return true
end

import Base.(+)
function +(i1::Impulse{T},i2::Impulse{T})::Impulse{T} where {T}
    in_time(t) = i1.in_time(t) + i2.in_time(t)
    in_freq(ω) = i1.in_freq(ω) + i2.in_freq(ω)

    Impulse{T}(in_time,in_freq)
end

import Base.(*)
function *(a,impulse::Impulse{T})::Impulse{T} where {T}

    if typeof(one(Complex{T})*a) != Complex{T}
        error("Multiplying impulse by $a would cause impulse type to change, please explicitly cast $a to same type as Source ($T)")
    end

    in_time(t) = a*impulse.in_time(t)
    in_freq(ω) = a*impulse.in_freq(ω)

    Impulse{T}(in_time,in_freq)
end

*(i::Impulse,a) = *(a,i)


"""
Delta function of unit area in the time domain, centred at t=t0.
"""
function TimeDeltaFunctionImpulse(t0::T) where {T<:AbstractFloat}
    in_time(t::T) = error("Time delta function is not defined in time")
    in_freq(ω::T) = exp(-im*ω*t0)
    Impulse{T}(in_time, in_freq, false)
end

"""
Delta function of unit area in the time domain, centred at t=t0.
"""
function FreqDeltaFunctionImpulse(ω0::T) where {T<:AbstractFloat}
    in_time(t::T) = exp(-im*ω0*t)
    in_freq(ω::T) = error("Frequency delta function is not defined in frequency")
    Impulse{T}(in_time, in_freq, false)
end

"""
Returns a gaussian impulse function in the time domain.
"""
function GaussianImpulse(maxω::T, a::T = T(2.48)/maxω^2) where {T<:AbstractFloat}
    in_time(t::T) = exp(-t^2/(4a))
    in_freq(ω::T) = exp(-a*ω^2)*(2sqrt(a*pi))
    Impulse{T}(in_time, in_freq)
end
