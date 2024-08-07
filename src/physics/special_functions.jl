export sbesselj, shankelh1, diffsbessel, diffbessel
export diffbesselj, diffhankelh1, diffsbesselj, diffshankelh1, diff3sbesselj, diff2sbesselj, diff2shankelh1, diff3shankelh1
export gaunt_coefficient
export associated_legendre_indices, spherical_harmonics_indices, lm_to_spherical_harmonic_index

export complex_legendre_array
export spherical_harmonics, spherical_harmonics_dθ
export cartesian_to_radial_coordinates, radial_to_cartesian_coordinates
export cartesian_to_spherical_coordinates, spherical_to_cartesian_coordinates
export spherical_to_cartesian_transform, cartesian_to_spherical_transform
export spherical_to_cartesian_vector, cartesian_to_spherical_vector
export atan

#NOTE should use https://github.com/JuliaMath/Bessels.jl instead of SpecialFunctions.jl

#NOTE atan(x,y) is defined here but will likely be added to base soon.

import Base.atan
atan(y::Complex,x::Complex) = - im * log( (x+y*im) / sqrt(x^2+y^2) )

"""
    sbesselj(m,x)

Returns the spherical besselj function. The order is 'm' and the argument is 'x'. Note 'x' can be a complex number.
"""
# function sbesselj(m::Number,x::T) where T<: Union{F,Complex{F}} where F <: AbstractFloat
    # if (abs(x) > eps(F))
    #     return sqrt(pi/(T(2)*x)) * besselj(m+1/2,x)
    # else
    #     return (m > 0 ? zero(T) : one(T))
    # end
# end
function sbesselj(nu::Number, x::T) where {T}
    besselj_nuhalf_x = besselj(nu + one(nu)/2, x)
    if abs(x) ≤ sqrt(eps(real(zero(besselj_nuhalf_x))))
        nu == 0 ? one(besselj_nuhalf_x) : zero(besselj_nuhalf_x)
    else
        √((float(T))(π)/2x) * besselj_nuhalf_x
    end
end

"""
    shankelh1(m,x)

Returns the spherical hankel function of the first kind. The order is 'm' and the argument is 'x'. Note 'x' can be a complex number.
"""
shankelh1(nu::Number,x::T) where T = √((float(T))(π)/2x) * hankelh1(nu+one(nu)/2,x)

"""
    diffsbessel(f::Function,m,x)

Differentiates the spherical bessel function 'f' (for any spherical bessel). The order is 'm' and the argument is 'x'. Note 'x' can be a complex number.
"""
function diffsbessel(f::Function,n::Number,z::Number)
    return f(n-1,z) - (n+1) * f(n,z) / z
end

# """
#     diffsbessel(f::Function,m,x,n::Int)

# Differentiates 'n' times any spherical Bessel function 'f' of order 'm' and at the argument 'x'. Uses the formula 
    
# ``(1/z d/dz)^n (z^{m+1} f_m(z)) = z^{m-n+1} f_{m-n}``

# which leads to     

# ``(1/z d/dz)^{n} (z^{m+1} f_m(z)) = (1/z d/dz)^{n-1} (1/z d/dz) (z^{m+1} f_m(z))``
    
# """
# function diffsbessel(f::Function,m::Number,z,n::Int)
#     if n == 0
#         return f(m, z)
#     elseif n > 0
#         n = n - 1
#         return diffsbessel(f,m,z)

        
#     else
#         error("Can not differentiate a negative number of times")
#     end
# end

function diffsbesselj(n::Number,z::Number)
    return if n == 0
        - sbesselj(1,z) # case due to numerical stability
    elseif z ≈ 0
        (n == 1) ? typeof(z)(1/3) : typeof(z)(0)
    else    
        sbesselj(n-1,z) - (n+1) * sbesselj(n,z) / z
    end
end

function diff2sbesselj(n::Number,z::Number)
    return if z ≈ 0
        if n == 0 
            -typeof(z)(1/3) 
        elseif n == 2
            -typeof(z)(2/15)
        else       
            typeof(z)(0)
        end
    else    
        2sbesselj(n+1,z) / z + (n^2 - n - z^2) * sbesselj(n,z) / z^2
    end
end

function diff2shankelh1(n::Number,z::Number)
    return 2shankelh1(n+1,z) / z + (n^2 - n - z^2) * shankelh1(n,z) / z^2
end

function diff3sbesselj(n::Number,z::Number)
    return ((-2 + n) * (-n + n^2 - z^2) * sbesselj(n,z) - 
    z * (6 + n + n^2 - z^2) * sbesselj(n+1,z))/z^3
end

function diff3shankelh1(n::Number,z::Number)
    return  ((-2 + n) * (-n + n^2 - z^2) * shankelh1(n,z) - 
    z * (6 + n + n^2 - z^2) * shankelh1(n+1,z))/z^3
end

function diffshankelh1(n::Number,z::Number)
    return shankelh1(n-1,z) - (n+1) * shankelh1(n,z) / z
end

"""
    diffbessel(f::Function,m,x,n::Int)

Differentiates 'n' times any bessel function 'f' of order 'm' and at the argument 'x'.
"""
function diffbessel(f::Function,m::Number,z,n::Int)
    if n == 0
        return f(m, z)
    elseif n > 0
        n = n - 1
        return 0.5*(diffbessel(f,m-1,z,n) - diffbessel(f,m+1,z,n))
    else
        error("Can not differentiate a negative number of times")
    end
end

"""
    diffhankelh1(m,x,n::Int)

Differentiates 'n' times the hankelh1 function of order 'm' and at the argument 'x'.
"""
function diffhankelh1(m::Number,z::T,n::Int) where T<:Number
    if n == 0
        return hankelh1(m, z)
    elseif n > 0
        n = n - 1
        return 0.5*(diffhankelh1(m-1,z,n) - diffhankelh1(m+1,z,n))
    else
        error("Can not differentiate a negative number of times")
    end
end

"Derivative of Hankel function of the first kind"
diffhankelh1(n,z) = 0.5*(hankelh1(-1 + n, z) - hankelh1(1 + n, z))

"Derivative of Bessel function of first kind"
diffbesselj(n::Number,z::T) where T<:Number = 0.5*(besselj(-1 + n, z) - besselj(1 + n, z))

"m-th Derivative of Hankel function of the first kind"
function diffbesselj(n::Number,z::T,m::Int) where T<:Number
    if m == 0
        return besselj(n, z)
    elseif m > 0
        m = m - 1
        return T(0.5)*(diffbesselj(n-1,z,m) - diffbesselj(n+1,z,m))
    else
        error("Can not differentiate a negative number of times")
    end
end


"""
    gaunt_coefficient(l1,m1,l2,m2,l3,m3)

A version of the Gaunt coefficients which are used to write the product of two spherical harmonics. If Y_{l,m} is a complex spherical harmonic, with the typical phase conventions from quantum mechanics, then:

    gaunt_coefficient(l1,m1,l2,m2,l3,m3) = 4*π*im^{l2+l3-l1} Integral[Y_{l1,m1}*conj(Y_{l2,m2})*conj(Y_{l3,m3})]

where the integral is over the solid angle.

The most standard gaunt coefficients `G(l1,m1;l2,m2;l3)` are related through the identity:

    4pi * G(l1,m1;l2,m2;l3) = im^(l1-l2-l3) * (-1)^m2 * gaunt_coefficient(l1,m1,l2,-m2,l3,m1+m2)

"""
function gaunt_coefficient(T::Type{<:AbstractFloat},l1::Int,m1::Int,l2::Int,m2::Int,l3::Int,m3::Int)
    # note the wigner3j has only one convention, and is highly symmetric.
    return (one(T)*im)^(l2+l3-l1) * (-T(1))^m1 * sqrt(4pi*(2*l1+1)*(2*l2+1)*(2*l3+1)) *
        wigner3j(T,l1,l2,l3,0,0,0) * wigner3j(T,l1,l2,l3,m1,-m2,-m3)
end
gaunt_coefficient(l1::Int,m1::Int,l2::Int,m2::Int,l3::Int,m3::Int) = gaunt_coefficient(Float64,l1,m1,l2,m2,l3,m3)

function lm_to_spherical_harmonic_index(l::Int,m::Int)::Int
    if l < abs(m)
        error("The order m of a spherical harmonic must be great or equal than the degree l.")
    end    
    return l^2 + m + l + 1
end    

function spherical_harmonics_indices(l_max::Int)
    ls = [l for l in 0:l_max for m in -l:l]
    ms = [m for l in 0:l_max for m in -l:l]

    return ls, ms
end

function associated_legendre_indices(l_max::Int)
    ls = [l for l in 0:l_max for m in 0:l]
    ms = [m for l in 0:l_max for m in 0:l]

    return ls, ms
end

"""
`complex_legendre_array(l_max::Int, z::Complex{T})`

returns a vector of all associated legendre functions with degree `l <= l_max`. The degree and order (indices) of the elements of the vector are given by `spherical_harmonics_indices(l_max::Int)`.

The associated legendre polynomials are taken from the package AssociatedLegendrePolynomials.jl.
"""
function complex_legendre_array(l_max::Int, z::Union{T,Complex{T}}) where {T<:AbstractFloat}

    leg_array = zeros(Complex{T}, (l_max + 1)^2)
    ind = 0
    for l in 0:l_max
        for m in -l:l
            ind += 1
            if m < 0
                leg_array[ind] = (-1)^m * Nlm(l, -m) * conj(Plm(l, -m, z))
            else
                leg_array[ind] = Nlm(l, m) * Plm(l, m, z)
            end
        end
    end

    return leg_array

end


"""
`spherical_harmonics(l_max::Int, θ::Complex{T}, φ::Union{T,Complex{T}})`

returns a vector of all spherical harmonics with degree `l <= l_max` with complex arguments. The degree and order (indices) of the elements of the vector are given by `spherical_harmonics_indices(l_max::Int)`.
"""
function spherical_harmonics(l_max::Int, θ::Complex{T}, φ::Union{T,Complex{T}}) where {T<:AbstractFloat}

    Plm_arr = complex_legendre_array(l_max, cos(θ))

    Ylm_vec = Vector{Complex{T}}(undef, (l_max + 1)^2)

    ind = 0
    for l in 0:l_max
        for m in -l:l
            ind += 1
            Ylm_vec[ind] = Plm_arr[ind] * exp(m * im * φ)
        end
    end

    return Ylm_vec
end


"""
    spherical_harmonics_dθ(l_max::Int, θ::T, φ::T)

Returns a vector of all ``\\partial Y_{(l,m)}(\\theta,\\phi) / \\partial \\theta`` for all the degrees `l` and orders `m`.
"""
function spherical_harmonics_dθ(l_max::Int, θ::T, φ::T) where T <: AbstractFloat

    c(l, m) = sqrt(Complex{T}(l^2 - m^2)) / sqrt(Complex{T}(4l^2 - 1))
    Ys = spherical_harmonics(l_max+1, θ, φ)

    lm_to_n = lm_to_spherical_harmonic_index

    dY1s = [
        (l == 0) ? zero(Complex{T}) : l * c(l + 1, m) * Ys[lm_to_n(l+1,m)] 
    for l = 0:l_max for m = -l:l]
    
    dY2s = [
        (l == abs(m)) ? zero(Complex{T}) : (l + 1) * c(l, m) * Ys[lm_to_n(l-1,m)]
    for l = 0:l_max for m = -l:l] 
    
    return (dY1s - dY2s) ./ sin(θ)
end    

"""
    spherical_harmonics(l_max::Int, θ::T, φ::T)

Returns a vector of all the spherial harmonics ``Y_{(l,m)}(\\theta,\\phi)`` for all the degrees `l` and orders `m`. The degree and order (indices) of the elements of the vector are given by `spherical_harmonics_indices(l_max::Int)`.
"""
function spherical_harmonics(l_max::Int, θ::T, φ::Union{T,Complex{T}}) where {T<:AbstractFloat}

    ls, ms = associated_legendre_indices(l_max)
    Plm_arr = sf_legendre_array(GSL_SF_LEGENDRE_SPHARM, l_max, cos(θ))[eachindex(ls)]

    Ylm_vec = Vector{Complex{T}}(undef, (l_max+1)^2)
    Ylm_vec[1] = Plm_arr[1]

    ind1 = 1
    ind2 = 1
    for i = 1:l_max
        inds1 = (ind1+i):(ind1+2i)
        Ylm_vec[(ind2+i):(ind2+3i)] = [reverse((-1).^ms[inds1[2:end]] .* conj(Plm_arr[inds1[2:end]])); Plm_arr[inds1]]
        ind1 += i
        ind2 += 2i
    end

    ls, ms = spherical_harmonics_indices(l_max)
    Ylm_vec = (-one(T)) .^ ms .* exp.(ms .* (im*φ)) .* Ylm_vec

    return Ylm_vec
end

cartesian_to_radial_coordinates(x::AbstractVector) = cartesian_to_radial_coordinates(SVector(x...))
radial_to_cartesian_coordinates(θ::AbstractVector) = radial_to_cartesian_coordinates(SVector(θ...))

cartesian_to_spherical_coordinates(x::AbstractVector) = cartesian_to_radial_coordinates(x)
spherical_to_cartesian_coordinates(θ::AbstractVector) = radial_to_cartesian_coordinates(θ)

spherical_to_cartesian_vector(v::AbstractVector,θ::AbstractVector) = spherical_to_cartesian_vector(SVector(v...),SVector(θ...))
cartesian_to_spherical_vector(v::AbstractVector,x::AbstractVector) = cartesian_to_spherical_vector(SVector(v...),SVector(x...))

spherical_to_cartesian_transform(θ::AbstractVector) = spherical_to_cartesian_transform(SVector(θ...))
cartesian_to_spherical_transform(x::AbstractVector) = cartesian_to_spherical_transform(SVector(x...))

function spherical_to_cartesian_transform(rθφ::SVector{3})
    r, θ, φ = rθφ
    M = [
        [sin(θ) * cos(φ)  cos(θ) * cos(φ) -sin(φ)];
        [sin(θ) * sin(φ)  cos(θ) * sin(φ)  cos(φ)];
        [cos(θ)  -sin(θ)  0]
    ]

    return M
end   

function cartesian_to_spherical_transform(xyz::SVector{3})
    r, θ, φ = cartesian_to_spherical_coordinates(xyz)
    M = [
        [ sin(θ) * cos(φ)  sin(θ) * sin(φ)  cos(θ)];
        [ cos(θ) * cos(φ)  cos(θ) * sin(φ)  -sin(θ)];
        [-sin(φ)  cos(φ)  0]
    ]

    return M
end   

function spherical_to_cartesian_vector(v::SVector{3}, rθφ::SVector{3})

    return spherical_to_cartesian_transform(rθφ) * v
end   

function cartesian_to_spherical_vector(v::SVector{3}, xyz::SVector{3}) 

    return cartesian_to_spherical_transform(xyz) * v
end   

function cartesian_to_radial_coordinates(x::SVector{3,CT}) where CT
    r = sqrt(sum(x .^2)) # note this is, and should be, a complex number when x is a complex vector
    θ = atan(sqrt(x[1]^2+x[2]^2),x[3])
    φ = (!iszero(x[1]) || !iszero(x[2])) ? atan(x[2], x[1]) : zero(CT)
    return [r,θ,φ]
end

function cartesian_to_radial_coordinates(x::SVector{2,CT}) where CT
    r = sqrt(sum(x .^2)) # note this should be complex if x is complex
    θ = atan(x[2], x[1])
    return [r,θ]
end

function radial_to_cartesian_coordinates(rθφ::SVector{3,CT}) where CT
    r, θ, φ = rθφ
    x = r * sin(θ) * cos(φ)
    y = r * sin(θ) * sin(φ)
    z = r * cos(θ)

    return [x,y,z]
end

function radial_to_cartesian_coordinates(rθ::SVector{2,CT}) where CT
    r, θ = rθ
    x = r * cos(θ)
    y = r * sin(θ)

    return [x,y]
end
