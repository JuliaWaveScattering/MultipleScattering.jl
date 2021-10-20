"""
    point_source(medium::Acoustic, source_position, amplitude=1)::Source

Create 2D [`Acoustic`](@ref) point [`Source`](@ref) (zeroth Hankel function of first type)
"""
function point_source(medium::Acoustic{T,2}, source_position::AbstractVector, amplitude::Union{T,Complex{T},Function} = one(T))::Source{T,Acoustic{T,2}} where T <: AbstractFloat

    # Convert to SVector for efficiency and consistency
    source_position = SVector{2,T}(source_position)

    if typeof(amplitude) <: Number
        amp(ω) = amplitude
    else
        amp = amplitude
    end
    source_field(x,ω) = (amp(ω)*im)/4 * hankelh1(0,ω/medium.c * norm(x-source_position))

    function source_coef(order,centre,ω)
        k = ω/medium.c
        r, θ = cartesian_to_radial_coordinates(centre - source_position)

        # using Graf's addition theorem
        return (amp(ω)*im)/4 * [hankelh1(-n,k*r) * exp(-im*n*θ) for n = -order:order]
    end

    return Source{T,Acoustic{T,2}}(medium, source_field, source_coef)
end

# If we replaced 3 with Dim below this could should work for all dimensions! Test carefully after changing.
function point_source(medium::Acoustic{T,3}, source_position, amplitude::Union{T,Complex{T},Function} = one(T))::Source{T,Acoustic{T,3}} where T <: AbstractFloat

    # Convert to SVector for efficiency and consistency
    source_position = SVector{3,T}(source_position)

    if typeof(amplitude) <: Number
        amp(ω) = amplitude
    else
        amp = amplitude
    end

    # source_field(x,ω) = amp(ω) / (T(4π) * norm(x-source_position)) * exp(im * ω/medium.c * norm(x-source_position))
    # source_field(x,ω) = amp(ω)/sqrt(4π) * shankelh1(0, ω/medium.c * norm(x-source_position))
    source_field(x,ω) = amp(ω) * outgoing_basis_function(medium, ω)(0,x-source_position)[1]

    # centre = 12.0 .* rand(3)
    # x = centre + rand(3)

    # order = 5
    # U = outgoing_translation_matrix(medium, order, ω,  centre);
    # vs = regular_basis_function(medium, ω)(order, x - centre);
    # us = outgoing_basis_function(medium, ω)(0, x);
    # sum(U[1,:] .* vs) - us[1]

    function source_coef(order,centre,ω)
        U = outgoing_translation_matrix(medium, order, ω,  centre);
        return amp(ω) * U[1,:]
    end

    return Source{T,Acoustic{T,3}}(medium, source_field, source_coef)
end


function plane_source(medium::Acoustic{T,Dim}; position::AbstractArray{T} = SVector(zeros(T,Dim)...),
        direction = SVector(one(T), zeros(T,Dim-1)...),
        amplitude::Union{T,Complex{T},Function} = one(T))::Source{T,Acoustic{T,Dim}} where {T, Dim}

    plane_source(medium, position, direction, amplitude)
end

"""
    plane_source(medium::Acoustic, source_position, source_direction=[1,0], amplitude=1)::Source

Create an [`Acoustic`](@ref) planar wave [`Source`](@ref)
"""
function plane_source(medium::Acoustic{T,2}, position::AbstractArray{T}, direction::AbstractArray{T} = SVector(one(T),zero(T)), amplitude::Union{T,Complex{T}} = one(T))::Source{T,Acoustic{T,2}} where {T}

    # Convert to SVector for efficiency and consistency
    position = SVector(position...)
    direction = SVector((direction ./ norm(direction))...) # unit direction

    # This pseudo-constructor is rarely called, so do some checks and conversions
    if iszero(norm(direction))
        throw(DomainError("Source direction must not have zero magnitude."))
    end

    if typeof(amplitude) <: Number
        amp(ω) = amplitude
    else
        amp = amplitude
    end

    source_field(x,ω) = amp(ω)*exp(im*ω/medium.c*dot(x-position, direction))

    function source_coef(order,centre,ω)
        # Jacobi-Anger expansion
        θ = atan(direction[2],direction[1])
        source_field(centre,ω) * [exp(im * n *(T(pi)/2 -  θ)) for n = -order:order]
    end

    return Source{T,Acoustic{T,2}}(medium, source_field, source_coef)
end

function plane_source(medium::Acoustic{T,3}, position::AbstractArray{T}, direction::AbstractArray{T} = SVector(zero(T),zero(T),one(T)), amplitude::Union{T,Complex{T}} = one(T)) where {T}

    # Convert to SVector for efficiency and consistency
    position = SVector(position...)
    direction = SVector( (direction ./ norm(direction))...) # unit direction

    # This pseudo-constructor is rarely called, so do some checks and conversions
    if iszero(norm(direction))
        throw(DomainError("Source direction must not have zero magnitude."))
    end

    if typeof(amplitude) <: Number
        amp(ω) = amplitude
    else
        amp = amplitude
    end

    source_field(x,ω) = amp(ω) * exp(im * ω/medium.c * dot(x-position, direction))

    function source_coef(order,centre,ω)
        # plane-wave expansion for complex vectors
        r, θ, φ  = cartesian_to_radial_coordinates(direction)
        Ys = spherical_harmonics(order, θ, φ)
        lm_to_n = lm_to_spherical_harmonic_index

        return T(4pi) * source_field(centre,ω) .*
        [
            Complex{T}(im)^l * (-one(T))^m * Ys[lm_to_n(l,-m)]
        for l = 0:order for m = -l:l]
    end

    return Source{T,Acoustic{T,3}}(medium, source_field, source_coef)
end

function regular_spherical_coefficients(psource::PlaneSource{T,Dim,1,Acoustic{T,Dim}}) where {T,Dim}
    source = plane_source(psource.medium;
        amplitude = psource.amplitude[1],
        position = psource.position,
        direction = psource.direction
    )

    source.coefficients
end
