
"""
Convert a FrequencySimulationResult into a TimeSimulationResult by using the inverse fourier transform.
Assumes only positive frequencies and a real time signal
"""
function frequency_to_time(simres::FrequencySimulationResult{T,Dim,FieldDim};
        t_vec = ω_to_t(simres.ω),
        impulse = GaussianImpulse(maximum(simres.ω)),
        method = :dft
    ) where {Dim,FieldDim,T}

    time_field = frequency_to_time(transpose(field(simres)), transpose(simres.ω), t_vec;
        impulse = impulse, method = method)

    return TimeSimulationResult(transpose(time_field), deepcopy(simres.x), t_vec)
end

"""
Convert a TimeSimulationResult into a FrequencySimulationResult by using the fourier transform.
Assumes only positive frequencies and a real time signal
"""
function time_to_frequency(timres::TimeSimulationResult{T,Dim,FieldDim};
        ω_vec = t_to_ω(timres.t),
        impulse = GaussianImpulse(maximum(ω_vec)),
        method =:dft
    ) where {Dim,FieldDim,T}

    freq_field = time_to_frequency(transpose(field(timres)), transpose(timres.t), ω_vec;
        impulse = impulse, method = method)

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
    t_arr = linspace(zero(T),2π/dω,2N+2)[1:(2N+1)]
    return t_arr
end

"The inverse of ω_to_t if ω_vec[1] == 0"
function t_to_ω(t_arr::AbstractArray{T}) where T <: AbstractFloat
    N = Int(round((length(t_arr)-one(T))/T(2)))
    maxt = t_arr[2] - t_arr[1] + t_arr[end]
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
Calculates the time response from the frequency response by approximating an
inverse Fourier transform. The time signal is assumed to be real and the
frequenices ω_vec are assumed to be positive (can include zero) and sorted. The
result is convoluted ωith the user specified impulse, which is a function of the
frequency.
"""
function frequency_to_time(field_mat::AbstractArray{Complex{T}}, ω_vec::AbstractVector{T},
        t_vec::AbstractArray{T} = ω_to_t(ω_vec);
        impulse::Impulse{T} = TimeDeltaFunctionImpulse(zero(T)),
        addzerofrequency=true, method=:dft) where T <: AbstractFloat
    # we assume the convention: f(t) = 1/(2π) ∫f(ω)exp(-im*ω*t)dω

    if size(field_mat,1) != size(ω_vec,1) error("Vector of frequencies ω_vec expected to be same size as size(field_mat,1)") end

    # approximate  ∫f(ω)exp(-im*ω*t)dω. The advantage of not using FFT is the ability to easily sample any time.
    impulse_in_freq_vec = map(impulse.in_freq, ω_vec)
    function f(t::T,j::Int)
        fs = impulse_in_freq_vec.*field_mat[:,j].*exp.(-(im*t).*ω_vec)
        if method == :dft && ω_vec[1] == zero(T)
            fs[1] = fs[1]/T(2) # so as to not add ω=0 tωice in the integral of ω over [-inf,inf]
        end
        fs
    end
    inverse_fourier_integral = (t,j) -> numerical_integral(ω_vec, f(t,j); method=method)
    u = [inverse_fourier_integral(t,j) for t in t_vec, j in indices(field_mat,2)]

    return real.(u)/pi # constant due to our Fourier convention and because only positive frequencies are used.
end

"""
The inverse of the function frequency_to_time (only an exact inverse when using
:dft integration)
"""
function time_to_frequency(field_mat::AbstractMatrix{T}, t_vec::AbstractVector{T},
        ω_vec::AbstractArray{T} = t_to_ω(t_vec);
        impulse::Impulse{T} = FreqDeltaFunctionImpulse(zero(T)),
        method=:dft) where T <: AbstractFloat

    # we assume the convention: f(ω) =  ∫f(t)exp(im*ω*t)dt
    impulse_in_time_vec = map(impulse.in_time, t_vec)
    f(ω::T, j::Int) = impulse_in_time_vec.*field_mat[:,j].*exp.((im*ω).*t_vec)
    fourier_integral = (ω,j) -> numerical_integral(t_vec, f(ω,j); method=method)
    uhat = [fourier_integral(ω,j) for ω in ω_vec, j in indices(field_mat,2)]

    return uhat
end

function numerical_integral(xs::AbstractArray{T}, fs::Union{AbstractArray{T},AbstractArray{Complex{T}}};
        method = :dft) where T <: AbstractFloat

    if method == :trapezoidal
        sum(1:(length(xs)-1)) do xi
            (fs[xi]+fs[xi+1])*(xs[xi+1] - xs[xi])/T(2)
        end
    elseif method == :dft
        fs[1]*(xs[2]-xs[1]) +
        sum(2:length(xs)) do xi
            fs[xi]*(xs[xi] - xs[xi-1])
        end
    else
        error("The method $method for numerical integration is not known.")
    end
end
