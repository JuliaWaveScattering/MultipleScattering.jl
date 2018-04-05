include("domain/particle.jl")

"Immutable struct"
struct Source
    field_expr::Expr
    coef_expr::Expr
end

import Base.(+)

function +(a::Source, b::Source)
    Source(
        parse(string(string(a.field_expr)," + ",string(b.field_expr))),
        parse(string(string(a.coef_expr)," + ",string(b.coef_expr))),
    )
end

"Helper function to return a plane wave source"
function planar_wave_source{T}(magnitude::Complex{T}, direction::Vector{T})
    # Have to use string interpolation
    Source(
        :((x,k)->($(magnitude))*exp(k*(x[1]*$(direction[1])+ x[2]*$(direction[2])))),
        :((p,k)->[1.0+0.0im]),
    )
end

"Safe access function for getting values of function from Source"
function get_source_field{T}(source::Source, x::Vector{T}, k::T)::Complex{T}
    field = eval(source.field_expr)
    field(x,k)
end

"Safe access function for getting coefficients from Source"
function get_coefs{T}(source::Source, p::Particle{T}, k::T)::Vector{Complex{T}}
    coefs = eval(source.coef_expr)
    coefs(p,k)
end
