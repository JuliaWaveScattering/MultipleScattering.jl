
"""
Convert a FrequencySimulationResult into a TimeSimulationResult by using the inverse fourier transform.
Assumes only positive frequencies and a real time signal
"""
function frequency_to_time(simres::FrequencySimulationResult{T,Dim,FieldDim};
        t_vec::AbstractVector{T} = ω_to_t(simres.ω),
        impulse::ContinuousImpulse{T} = TimeDiracImpulse(zero(T)), #GaussianImpulse(maximum(simres.ω)),
        discrete_impulse::DiscreteImpulse{T} = continuous_to_discrete_impulse(impulse, t_vec, simres.ω),
        method = :dft
    ) where {Dim,FieldDim,T}

    time_field = frequency_to_time(transpose(field(simres)), simres.ω, t_vec;
        discrete_impulse = discrete_impulse, method = method)

    return TimeSimulationResult(transpose(time_field), simres.x, t_vec)
end

"""
Convert a TimeSimulationResult into a FrequencySimulationResult by using the fourier transform.
Assumes only positive frequencies and a real time signal
"""
function time_to_frequency(timres::TimeSimulationResult{T,Dim,FieldDim};
        ω_vec = t_to_ω(timres.t),
        impulse::ContinuousImpulse{T} = TimeDiracImpulse(zero(T)), #GaussianImpulse(maximum(ω_vec)),
        discrete_impulse::DiscreteImpulse{T} = continuous_to_discrete_impulse(impulse,timres.t, ω_vec),
        method =:dft
    ) where {Dim,FieldDim,T}

    freq_field = time_to_frequency(transpose(field(timres)), timres.t, ω_vec;
        discrete_impulse = discrete_impulse, method = method)

    return FrequencySimulationResult(transpose(freq_field), deepcopy(timres.x), ω_vec)
end

"""
returns an array of time from the frequency array ω_vec.
Uses the same convention for sampling the time as the discrete Fourier transfrom.
Assumes ω_vec is ordered and non-negative.
"""
function ω_to_t(ω_arr::AbstractArray{T}) where T <: AbstractFloat
    N = length(ω_arr)
    if ω_arr[1] == zero(T)
        N -= 1
    elseif minimum(ω_arr) < zero(T)
        error("expected only non-negative values for the frequencies")
    end
    dω = ω_arr[2] - ω_arr[1]
    t_arr = LinRange(zero(T),2π/dω,2N+2)[1:(2N+1)]
    return t_arr
end

"The inverse of ω_to_t if ω_vec[1] == 0"
function t_to_ω(t_arr::AbstractArray{T}) where T <: AbstractFloat
    N = Int(round((length(t_arr)-one(T))/T(2)))
    maxt = t_arr[2] - t_arr[1] + t_arr[end] - t_arr[1] # subtract t_arr[1] in case t_arr[1] != zero(T)
    maxω = N*2π/maxt
    ω_vec = LinRange(zero(T),maxω,N+1)
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
See also: [`DiscreteImpulse`](@ref), [`ContinuousImpulse`](@ref)

Calculates the time response from the frequency response by approximating an
inverse Fourier transform. The time signal is assumed to be real and the
frequenices ω_vec are assumed to be positive (can include zero) and sorted. The
result is convoluted in time ωith the user specified impulse.

We use the Fourier transform convention:
F(ω) =  ∫ f(t)*exp(im*ω*t) dt
f(t) = (2π)^(-1) * ∫ F(ω)*exp(-im*ω*t) dt

To easily sample any time, the default is not FFT, but a discrete version of the transform above.
"""
function frequency_to_time(field_mat::AbstractArray{Complex{T}}, ω_vec::AbstractVector{T},
        t_vec::AbstractArray{T} = ω_to_t(ω_vec);
        impulse::ContinuousImpulse{T} = TimeDiracImpulse(zero(T)),
        discrete_impulse::DiscreteImpulse{T} = continuous_to_discrete_impulse(impulse, t_vec, ω_vec),
        method=:dft) where T <: AbstractFloat

    if size(field_mat,1) != size(ω_vec,1) error("Vector of frequencies ω_vec expected to be same size as size(field_mat,1)") end

    function f(t::T,j::Int)
        fs = discrete_impulse.in_freq.*field_mat[:,j].*exp.(-(im*t).*ω_vec)
        if method == :dft && ω_vec[1] == zero(T)
            fs[1] = fs[1]/T(2) # so as to not add ω=0 tωice in the integral of ω over [-inf,inf]
        end
        fs
    end
    inverse_fourier_integral = (t,j) -> numerical_integral(ω_vec, f(t,j), method)
    u = [inverse_fourier_integral(t,j) for t in t_vec, j in axes(field_mat,2)]

    return real.(u)/pi # a constant 1/(2pi) appears due to our Fourier convention, but because we only use positive frequencies, and assume a real time signal, this becomes 1/pi.
end

"""
The inverse of the function frequency_to_time (only an exact inverse when using
:dft integration). We use the Fourier transform convention:
F(ω) =  ∫ f(t)*exp(im*ω*t) dt
"""
function time_to_frequency(field_mat::AbstractArray{T}, t_vec::AbstractVector{T},
        ω_vec::AbstractArray{T} = t_to_ω(t_vec);
        impulse::ContinuousImpulse{T} = TimeDiracImpulse(zero(T)),
        discrete_impulse::DiscreteImpulse{T} = continuous_to_discrete_impulse(impulse, t_vec, ω_vec),
        method=:dft) where T <: AbstractFloat

    # to use an impulse below in time we would need to do a discrete convolution, which we decided against.
    f(ω::T, j::Int) = field_mat[:,j].*exp.((im*ω).*t_vec)
    fourier_integral = (ω,j) -> numerical_integral(t_vec, f(ω,j), method)
    uhat = [discrete_impulse.in_freq[i]*fourier_integral(ω_vec[i],j) for i in eachindex(ω_vec), j in axes(field_mat,2)]

    return uhat
end

function numerical_integral(xs::AbstractArray{T}, fs::Union{AbstractArray{T},AbstractArray{Complex{T}}}, method = :dft) where T <: AbstractFloat

    if method == :trapezoidal
        sum(1:(length(xs)-1)) do xi
            (fs[xi] + fs[xi+1])*(xs[xi+1] - xs[xi])/T(2)
        end
    elseif method == :dft
        fs[1]*(xs[2] - xs[1]) +
        sum(2:length(xs)) do xi
            fs[xi]*(xs[xi] - xs[xi-1])
        end
    else
        error("The method $method for numerical integration is not known.")
    end
end
