
function besselj_field(source::Source{Acoustic{T,2},T}, medium::Acoustic{T,2}, centre::AbstractVector{T}; basis_order = 4) where T

    field(x,ω) = sum(
        source.coef(n,centre,ω)*besselj(n,ω/medium.c*norm(x - centre))*exp(im*n*atan2(x[2] - centre[2],x[1] - centre[1]))
    for n = -basis_order:basis_order)

   return field
end

function TwoDimAcousticPointSource{T}(medium::Acoustic{T,2}, source_position::AbstractVector{T}, amplitude::T = one(T))::Source{Acoustic{T,2},T}
    field(x::AbstractVector{T},ω) where T = amplitude*Complex{T}(im/4)*hankelh1(0,ω/medium.c*norm(x-source_position))
    # using Graf's addition theorem
    coef(n,centre,ω) = amplitude*Complex{T}(im/4)*hankelh1(-n,ω/medium.c*norm(centre - source_position))*
        exp(-Complex{T}(im)*n*atan2(centre[2] - source_position[2], centre[1] - source_position[1]))

    return Source{Acoustic{T,2},T}(field,coef)
end

function TwoDimAcousticPlanarSource{T}(medium::Acoustic{T,2}, source_position::AbstractVector{T}, source_direction::AbstractVector{T} = [one(T),zero(T)], amplitude::T = one(T))::Source{Acoustic{T,2},T}
    field(x::AbstractVector{T},ω) where T = amplitude*exp(im*ω/medium.c*dot(x-source_position,source_direction))
    # Jacobi-Anger expansion
    coef(n,centre,ω) = field(centre,ω) * exp(im * n *(T(pi)/2 - atan2(source_direction[2],source_direction[1]) ))

    return Source{Acoustic{T,2},T}(field,coef)
end
