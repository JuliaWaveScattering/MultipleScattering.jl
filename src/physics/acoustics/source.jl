
function besselj_field(source::Source{Acoustic{T,2},T}, medium::Acoustic{T,2}, centre; basis_order = 4) where T

    # Convert to SVector for efficiency and consistency
    centre = SVector{2,T}(centre)

    return (x,ω) -> sum(
        source.coef(n,centre,ω)*besselj(n,ω/medium.c*norm(x - centre))*exp(im*n*atan2(x[2] - centre[2],x[1] - centre[1]))
    for n = -basis_order:basis_order)

end

function point_source{T}(medium::Acoustic{T,2}, source_position, amplitude::Union{T,Complex{T}} = one(T))::Source{Acoustic{T,2},T}

    # Convert to SVector for efficiency and consistency
    source_position = SVector{2,T}(source_position)

    source_field(x,ω) = (amplitude*im)/4*hankelh1(0,ω/medium.c*norm(x-source_position))

    function source_coef(n,centre,ω)
        k = ω/medium.c
        r = norm(centre - source_position)
        θ = atan2(centre[2]-source_position[2], centre[1]-source_position[1])
        # using Graf's addition theorem
        return (amplitude*im)/4 * hankelh1(-n,k*r) * exp(-im*n*θ)
    end

    return Source{Acoustic{T,2},T}(source_field, source_coef)
end

function plane_source(medium::Acoustic{T,2}, source_position, source_direction = SVector(one(T),zero(T)), amplitude::Union{T,Complex{T}} = one(T))::Source{Acoustic{T,2},T} where {T}

    # Convert to SVector for efficiency and consistency
    source_position = SVector{2,T}(source_position)
    source_direction = SVector{2,T}(source_direction./norm(source_direction)) # unit direction

    # This pseudo-constructor is rarely called, so do some checks and conversions
    if iszero(norm(source_direction))
        throw(DomainError("Source direction must not have zero magnitude."))
    end

    source_field(x,ω) = amplitude*exp(im*ω/medium.c*dot(x-source_position, source_direction))

    function source_coef(n,centre,ω)
        # Jacobi-Anger expansion
        θ = atan2(source_direction[2],source_direction[1])
        source_field(centre,ω) * exp(im * n *(T(pi)/2 -  θ))
    end

    return Source{Acoustic{T,2},T}(source_field, source_coef)
end
