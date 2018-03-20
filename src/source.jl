
type Source{P,T} where P <: PhysicalProperties{T},T <: AbstractFloat
    get_field::function
    get_coefs::function
end

function AcousticPointSource{T}(source_position::Vector{T}, amplitude::T)::Source{Acoustics{T},T}
    get_field(x,k) = amplitude*hankel(k*norm(x-source_position))
    get_coefs(center,radius,k) = amplitude*hankel(k*norm(x-source_position))

    return Source{Acoustics{T},T}()
end
