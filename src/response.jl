## This is where the actual maths happens

function response{T}(model::FrequencyModel{T}, k::T)
    coefficient_matrix = scattering_coefficients_matrix(model, k)
    # the matrix of scattering coefficient mat is given by mat[m,p] = Zn_vec[m,p]*A_mat[m,p]
    # where mat[m+hankel_order+1,p] is the coefficient for hankelh1[m, kr] of the p-th particle[:,p]

    return [response_at(model, k, coefficient_matrix, model.listener_positions[:,i]) for i in 1:size(model.listener_positions,2)]
end

function response_at{T}(model::FrequencyModel{T}, k::T, coefficient_matrix::Matrix{Complex{T}}, x::Vector{T})
    # check if x is inside a scatterer
    inside_p(p) = norm(x - p.x) < p.r
    i = findfirst(inside_p, model.particles)
    # one(Complex{T}) is the incident wave when emitted and measured at x
    if i ==0
        return exp(im * k * dot(x - model.source_position, model.source_direction)) +
          scattered_waves_at(model, k, coefficient_matrix, x)
        # one(Complex{T}) + scattered_waves_at(model, k, coefficient_matrix, x) * exp(-im * k * dot(x, model.source_direction))
    else
        return internal_wave_at(model, k, coefficient_matrix, x, i)
    end
end

"a scattered waves is the total response minus the incident"
function scattered_waves_at{T}(model, k::T, coefficient_matrix::Matrix{Complex{T}}, x::Vector{T})
    Nh = model.hankel_order
    krs = [k*norm(x - p.x) for p in model.particles]
    θs  = [atan2(x[2] - p.x[2], x[1] - p.x[1]) for p in model.particles]

    scattered_from(i) = sum(-Nh:Nh) do m
          coefficient_matrix[m+Nh+1, i] * hankelh1(m, krs[i]) * exp(im * m * θs[i])
      end
    #Sum wave scattered from each particle
    return sum(scattered_from, eachindex(model.particles))
end

"the internal wave exists only within a particle and is the total response"
function internal_wave_at{T}(model::FrequencyModel{T}, k::T, coefficient_matrix::Matrix{Complex{T}}, x::Vector{T}, i::Int)
    # B(m,i) = A[m,i] Zn(model,p,k,m)/ Jm(k_j a_j)(Hm(k a_j) - Jm(k a_j)/Zn(model,p,k,m) )
    Nh = model.hankel_order
    p = model.particles[i]
    γ = model.c/p.c
    # If the particle is impenetrable then return zero(T)
    if abs(p.c) == T(Inf) || abs(p.ρ) == T(Inf) || abs(model.ρ) == zero(T)
      return zero(T)
    end
    internal_coefficient(m) = (coefficient_matrix[m+Nh+1,i]/besselj(m, γ*k*p.r))*
        (hankelh1(m,k*p.r) - besselj(m,k*p.r)/Zn(model,p,k,m))

    R = norm(x-p.x)
    Θ = atan2(x[2] - p.x[2], x[1] - p.x[1])
    return sum(-Nh:Nh) do m
            internal_coefficient(m)*besselj(m,γ*k*R)*exp(im*m*Θ)
      end

end

function scattering_coefficients_matrix{T}(model::FrequencyModel{T}, k::T)::Matrix{Complex{T}}
    # Generate response for one specific k
    # Number of particles
    P = length(model.particles)
    # No particles means wave is unchanged, response matrix is a scalar
    if P == 0
        return one(Complex{T})
    end
    # Highest hankel order used
    Nh = model.hankel_order
    # Number of hankel basis function at each particle
    H = 2Nh + 1
    # Total hankel basis functions in whole model
    B = H * P

    # hankel_order is the order of Hankel function (2hankel_order +1 = # of unknowns) for each particle
    # deciding the type and strength of the scatterer is in the function Zn.
    Zn_mat = [Zn(model,p,k, m) for m=-model.hankel_order:model.hankel_order, p in model.particles]

    # Pre-calculate these to save doing it for all n and m
    # matrix of angles of the vectors pointing from particle p to particle s
    θ_mat = [atan2(s.x[2] - p.x[2], s.x[1] - p.x[1]) for p in model.particles, s in model.particles]

    # k times distance from particle p to particle s
    kr_mat = [k*norm(s.x .- p.x) for p in model.particles, s in model.particles]

    # Function we will use to generate the L matrix
    function L_fnc(m::Int, n::Int, s::Int, p::Int)
        if s == p
            zero(Complex{T})
        else
            hankelh1(T(n-m), kr_mat[p, s]) * Zn_mat[n+Nh+1,p] * exp(im * (n - m) * θ_mat[p, s])
        end
    end
    # Generate \mathbb{L} or LL
    L_mat = [L_fnc(m, n, s, p) for m=-Nh:Nh, s=1:P, n=-Nh:Nh, p=1:P]
    L_mat = reshape(L_mat, (B, B))

    # Vector which represents the inncident wave in the k direction i.e. exp(im* dot(k,x)) for each basis function
    ψ = reshape([exp(im * (k * dot(p.x - model.source_position, model.source_direction) + m * T(pi) / 2)) for m=-Nh:Nh, p in model.particles], B)
    # the scattering coefficients A_mat in matrix, A[1,2] multiplies hankelh1(1, k*r) contributing to the 2nd particle scattered wave
    A_mat = reshape(-(L_mat + eye(T,B)) \ ψ, (H, P))

    return Zn_mat .* A_mat
end

"Returns a ratio used in multiple scattering which reflects the material properties of the particles"
function Zn{T}(model::FrequencyModel{T}, p::Particle{T}, k::T,  m::Int)
    m = T(abs(m))
    ak = p.r*k
    # check for material properties that don't make sense or haven't been implemented
    if abs(p.c*p.ρ) == T(NaN)
        error("scattering from a particle with density =$(p.ρ) and phase speed =$(p.c) is not defined")
    elseif abs(model.c*model.ρ) == T(NaN)
        error("wave propagation in a medium with density =$(model.ρ) and phase speed =$(model.c) is not defined")
    elseif abs(model.c) == zero(T)
        error("wave propagation in a medium with phase speed =$(model.c) is not defined")
    elseif abs(model.ρ) == zero(T) && abs(p.c*p.ρ) == zero(T)
        error("scattering in a medium with density $(model.ρ) and a particle with density =$(p.ρ) and phase speed =$(p.c) is not defined")
    end

    # set the scattering strength and type
    if abs(p.c) == T(Inf) || abs(p.ρ) == T(Inf)
        numer = diffbesselj(m, ak)
        denom = diffhankelh1(m, ak)
    elseif abs(model.ρ) == zero(T)
        γ = model.c/p.c #speed ratio
        numer = diffbesselj(m, ak) * besselj(m, γ * ak)
        denom = diffhankelh1(m, ak) * besselj(m, γ * ak)
    else
        q = (p.c*p.ρ)/(model.c*model.ρ) #the impedance
        γ = model.c/p.c #speed ratio
        numer = q * diffbesselj(m, ak) * besselj(m, γ * ak) - besselj(m, ak)*diffbesselj(m, γ * ak)
        denom = q * diffhankelh1(m, ak) * besselj(m, γ * ak) - hankelh1(m, ak)*diffbesselj(m, γ * ak)
    end

    return numer / denom
end

"Derivative of Hankel function of zeroth kind"
function diffhankelh1(n, z)
    if n != 0
        (hankelh1(-1 + n, z) - hankelh1(1 + n, z)) / 2
    else
        -hankelh1(typeof(n)(1), z)
    end
end

"Derivative of Bessel function of first kind"
function diffbesselj(n, z)
    if n != 0
        (besselj(-1 + n, z) - besselj(1 + n, z)) / 2
    else
        -besselj(typeof(n)(1), z)
    end
end
