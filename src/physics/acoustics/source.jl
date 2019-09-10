
function besselj_field(source::Source{Acoustic{T,2},T}, medium::Acoustic{T,2}, centre; basis_order = 4) where T

    # Convert to SVector for efficiency and consistency
    centre = SVector{2,T}(centre)

    return (x,ω) -> sum(
        source.coef(n,centre,ω)*besselj(n,ω/medium.c*norm(x - centre))*exp(im*n*atan(x[2] - centre[2],x[1] - centre[1]))
    for n = -basis_order:basis_order)

end

"""
    point_source(medium::Acoustic, source_position, amplitude=1)::Source{Acoustic}

Create 2D [`Acoustic`](@ref) point [`Source`](@ref) (zeroth Hankel function of first type)
"""
function point_source(medium::Acoustic{T,2}, source_position, amplitude::Union{T,Complex{T},Function} = one(T))::Source{Acoustic{T,2},T} where T <: AbstractFloat

    # Convert to SVector for efficiency and consistency
    source_position = SVector{2,T}(source_position)

    if typeof(amplitude) <: Number
        amp(ω) = amplitude
    else
        amp = amplitude
    end
    source_field(x,ω) = (amp(ω)*im)/4*hankelh1(0,ω/medium.c*norm(x-source_position))

    function source_coef(n,centre,ω)
        k = ω/medium.c
        r = norm(centre - source_position)
        θ = atan(centre[2]-source_position[2], centre[1]-source_position[1])
        # using Graf's addition theorem
        return (amp(ω)*im)/4 * hankelh1(-n,k*r) * exp(-im*n*θ)
    end

    return Source{Acoustic{T,2},T}(medium, source_field, source_coef)
end

function plane_source(medium::Acoustic{T,2}; position = SVector(zero(T),zero(T)),
        direction = SVector(one(T),zero(T)),
        amplitude::Union{T,Complex{T},Function} = one(T))::Source{Acoustic{T,2},T} where {T}

    plane_source(medium, position, direction, amplitude)
end

"""
    plane_source(medium::Acoustic, source_position, source_direction=[1,0], amplitude=1)::Source{Acoustic}

Create 2D [`Acoustic`](@ref) planar wave [`Source`](@ref)
"""
function plane_source(medium::Acoustic{T,2}, position, direction = SVector(one(T),zero(T)),
        amplitude::Union{T,Complex{T}} = one(T))::Source{Acoustic{T,2},T} where {T}

    # Convert to SVector for efficiency and consistency
    position = SVector{2,T}(position)
    direction = SVector{2,T}(direction./norm(direction)) # unit direction

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

    function source_coef(n,centre,ω)
        # Jacobi-Anger expansion
        θ = atan(direction[2],direction[1])
        source_field(centre,ω) * exp(im * n *(T(pi)/2 -  θ))
    end

    return Source{Acoustic{T,2},T}(medium, source_field, source_coef)
end
