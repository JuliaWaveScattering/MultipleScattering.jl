# This is where the actual maths happens

function response{T}(model::FrequencyModel{T}, k::T)
    
    response_matrix = build_response_matrix(model, k)
    # the matrix of scattering coefficient mat is given by mat[m,p] = Zn_vec[m]*A_mat[m,p]
    # where mat[m+hankel_order+1,p] is the coefficient for hankelh1[m, kr] of the p-th particle[:,p]
    #(the 1.0+0.0im is the incident wave when emitted and measured at listener_positions)
    return one(Complex{T}) + scattered_wave(model, k, response_matrix)
end

function build_response_matrix{T}(model::FrequencyModel{T}, k::T)
    # Generate response for one specific k
    # Number of particles
    P = length(model.particles)
    # No particles means wave is unchanged, response matrix is a scalar
    if P == 0
        return one(Complex{T})
    end
    
    a = model.particles[1].r # assume all have the same radius
    # Highest hankel order used
    Nh = model.hankel_order
    # Number of hankel basis function at each particle
    H = 2Nh + 1
    # Total hankel basis functions in whole model
    B = H * P
    
    # hankel_order is the order of Hankel function (2hankel_order +1 = # of unknowns) for each particle and using Zn(k*a,m,0,1.) gives Dirichlet boundary conditions
    Zn_vec = [ Zn(k*a, m, 0, 1.0) for m=-model.hankel_order:model.hankel_order]

    # Pre-calculate these to save doing it for all n and m
    # matrix of angles of the vectors pointing from particle p to particle s
    θ_mat = [atan2(s.x[2] - p.x[2], s.x[1] - p.x[1]) for p in model.particles, s in model.particles]
    # println("Built θ_mat, $(summary(θ_mat))")
    # k times distance from particle p to particle s
    kr_mat = [k*norm(s.x .- p.x) for p in model.particles, s in model.particles]
    # println("Built kr_mat, $(summary(kr_mat))")
    # Function we will use to generate the L matrix
    function L_fnc(m::Int, n::Int, s::Int, p::Int)
        if s == p
            zero(Complex{T})
        else
            Complex{T}(hankelh1(n-m, kr_mat[p, s])) * Zn_vec[n+Nh+1] * exp(im * (n - m) * θ_mat[p, s])
        end
    end
    # Generate \mathbb{L} or LL
    L_mat = [L_fnc(m, n, s, p) for m=-Nh:Nh, s=1:P, n=-Nh:Nh, p=1:P]
    # println("Built L_mat, $(summary(L_mat))")
    L_mat = reshape(L_mat, (B, B))
    
    # Vector which represents the ncident wave in k direction i.e. exp(im* dot(k,x)) for each basis function
    ψ = reshape([exp(im * (k * dot(p.x, model.source_direction) + m * pi / 2)) for m=-Nh:Nh, p in model.particles], B)


    # the scattering coefficients A_mat in matrix, A[1,2] multiplies hankelh1(1, k*r) contributing to the 2nd particle scattered wave
    A_mat = reshape(-(L_mat + eye(B)) \ ψ, (H, P))

    return Zn_vec .* A_mat
end

"a scattered wave is the total response minus the incident"
function scattered_wave{T}(model::FrequencyModel{T}, k::T, response_matrix::Matrix{Complex{T}})
    # Number of particles
    P = length(model.particles)
    # assume all particles have the same radius
    a = model.particles[1].r
    # Highest hankel order used
    Nh = model.hankel_order
    # Number of hankel basis function at each particle
    H = 2Nh + 1
    # Total hankel basis functions in whole model
    B = H * P
    
    #scattered wave from particles[p] measured at x
    function outwave(i::Int, x::Vector{T})
        p = model.particles[i]
        kr = k*norm(x - p.x)
        θ = atan2(x[2] - p.x[2], x[1] - p.x[1])
        resp_m(m) = response_matrix[m+Nh+1, i] * hankelh1(m, kr) * exp(im * m * θ)
        return sum(resp_m, -Nh:Nh)
    end

    function scattered_at(x)
        # Check it isn't inside a particle
        for p in model.particles
            if(norm(x - p.x) + 10 * eps(T) < a )
                warn("Listener at $x is inside the scatterer centred at $(p.x), this has not been implemented")
                return -one(Complex{T})
            end
        end
        # If measuring point is within particle then return zero else divide all the scattered waves measured at X by the incident wave
        outwave_at_x(i) = outwave(i, x)
        return sum(outwave_at_x,1:P) * exp(-im * k * dot(x, model.source_direction))
    end

    [scattered_at(model.listener_positions[:,i]) for i in 1:size(model.listener_positions,2)]
end

"Returns a ratio used in multiple scattering which reflects the material properties of the particles"
function Zn(ak::Float64, m::Int, q=Inf, γ=1.0)
    m = abs(m)
    if q != Inf && γ === Inf
        error("Haven't done the maths yet for γ = Inf and q != Inf")
    end

    if q === Inf
        numer = diffbesselj(m, ak)
        denom = diffhankelh1(m, ak)
    else
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
        -hankelh1(1, z)
    end
end

"Derivative of Bessel function of first kind"
function diffbesselj(n, z)
    if n != 0
        (besselj(-1 + n, z) - besselj(1 + n, z)) / 2
    else
        -besselj(1, z)
    end
end
