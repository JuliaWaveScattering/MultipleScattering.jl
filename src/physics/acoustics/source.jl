"""
    point_source(medium::Acoustic, source_position, amplitude=1)::RegularSource

Create 2D [`Acoustic`](@ref) point [`RegularSource`](@ref) (zeroth Hankel function of first type)
"""
function point_source(medium::Acoustic{T,2}, source_position::AbstractVector, amplitude::Union{T,Complex{T},Function} = one(T))::RegularSource{Acoustic{T,2}} where T <: AbstractFloat

    # Convert to SVector for efficiency and consistency
    source_position = SVector{2,T}(source_position)

    amp = amplitude isa Function ? amplitude : ω -> amplitude

    source_field(x,ω) = (amp(ω)*im)/4 * hankelh1(0,ω/medium.c * norm(x-source_position))

    function source_coef(order,centre,ω)
        k = ω/medium.c
        r, θ = cartesian_to_radial_coordinates(centre - source_position)

        # using Graf's addition theorem
        return (amp(ω)*im)/4 * [hankelh1(-n,k*r) * exp(-im*n*θ) for n = -order:order]
    end

    return RegularSource{Acoustic{T,2},WithoutSymmetry{2}}(medium, source_field, source_coef)
end

# If we replaced 3 with Dim below this could should work for all dimensions! Test carefully after changing.
function point_source(medium::Acoustic{T,3}, source_position, amplitude::Union{T,Complex{T},Function} = one(T))::RegularSource{Acoustic{T,3}} where T <: AbstractFloat

    # Convert to SVector for efficiency and consistency
    source_position = SVector{3,T}(source_position)

    amp = amplitude isa Function ? amplitude : ω -> amplitude

    source_field(x,ω) = amp(ω) * outgoing_basis_function(medium, ω)(0,x-source_position)[1]

    function source_coef(order,centre,ω)
        # translate the outgoing wave of order 0 to the centre
        U = outgoing_translation_matrix(medium, order, 0, ω, centre - source_position)
        return amp(ω) .* U[1,:]
    end

    return RegularSource{Acoustic{T,3},WithoutSymmetry{3}}(medium, source_field, source_coef)
end


function plane_source(medium::Acoustic{T,Dim}; position::AbstractArray{T} = SVector(zeros(T,Dim)...),
        direction = SVector(one(T), zeros(T,Dim-1)...),
        amplitude::Union{T,Complex{T},Function} = one(T),
        kws...
    )::RegularSource{Acoustic{T,Dim}} where {T, Dim}

    plane_source(medium, position, direction, amplitude; kws...)
end

"""
    plane_source(medium::Acoustic, position, direction = azimuthalnormal(Dim), amplitude = 1; causal = false, beam_width = Inf)::RegularSource

Create an [`Acoustic`](@ref) planar wave [`RegularSource`](@ref). The default `direction` is `[1,0]` in 2D and `[0,0,1]` in 3D. The `amplitude` can be a number or a function of the angular frequency.
"""
function plane_source(medium::Acoustic{T,Dim}, position::AbstractArray{T},
        direction::AbstractArray{T} = SVector{Dim,T}(azimuthalnormal(Dim)),
        amplitude::Union{T,Complex{T},Function} = one(T);
        causal::Bool = false,
        beam_width::T = T(Inf)
    )::RegularSource{Acoustic{T,Dim}} where {T,Dim}

    if iszero(norm(direction))
        throw(DomainError("RegularSource direction must not have zero magnitude."))
    end

    # Convert to SVector for efficiency and consistency
    position = SVector{Dim,T}(position)
    direction = SVector{Dim,T}(direction ./ norm(direction)) # unit direction

    S = (abs(dot(direction,azimuthalnormal(Dim))) == one(T)) ? PlanarAzimuthalSymmetry{Dim} : PlanarSymmetry{Dim}

    amp = amplitude isa Function ? amplitude : ω -> amplitude

    function source_field(x,ω)
        x_width = norm((x - position) - dot(x - position, direction)*direction)
        if (causal && dot(x - position,direction) < 0) || x_width > beam_width/2
            zero(Complex{T})
        else
            amp(ω)*exp(im*ω/medium.c * dot(x-position, direction))
        end
    end

    source_coef(order,centre,ω) = source_field(centre,ω) .* plane_wave_coefficients(direction, order)

    return RegularSource{Acoustic{T,Dim},S}(medium, source_field, source_coef)
end

# The coefficients of a plane wave with unit `direction` expanded in the regular basis
function plane_wave_coefficients(direction::SVector{2,T}, order::Integer) where T
    # Jacobi-Anger expansion
    θ = atan(direction[2],direction[1])
    return [exp(im * n *(T(pi)/2 -  θ)) for n = -order:order]
end

function plane_wave_coefficients(direction::SVector{3,T}, order::Integer) where T
    # plane-wave expansion for complex vectors
    r, θ, φ  = cartesian_to_radial_coordinates(direction)
    Ys = spherical_harmonics(order, θ, φ)
    lm_to_n = lm_to_spherical_harmonic_index

    return [
        T(4pi) * Complex{T}(im)^l * (-one(T))^m * Ys[lm_to_n(l,-m)]
    for l = 0:order for m = -l:l]
end

function regular_spherical_coefficients(psource::PlaneSource{T,Dim,1,Acoustic{T,Dim}}) where {Dim,T}

    source = plane_source(psource.medium;
        amplitude = psource.amplitude[1],
        position = psource.position,
        direction = psource.direction
    )

    source.coefficients
end
