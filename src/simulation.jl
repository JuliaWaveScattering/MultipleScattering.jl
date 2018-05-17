abstract type Simulation{Dim,T} end

mutable struct FrequencySimulation{Dim,P<:PhysicalProperties,T<:AbstractFloat} <: Simulation{Dim,T}
    medium::P
    particles::Vector{AbstractParticle{Dim,T}}
    source::Source{P,T}
end


import Base.run

# Main run function, all other run functions use this
function run(sim::FrequencySimulation{Dim,P,T}, ω::T, x_vec::Vector{SVector{Dim,T}}; basis_order::Int = 5) where {Dim,FieldDim,T,P<:PhysicalProperties{Dim,FieldDim,T}}

    # Calculate the Hankel coefficients around each particle, this is where most of the maths happens
    a_vec = basis_coefficients(sim, ω; basis_order=basis_order)

    # Evaluate the total field at the requested x positions
    field_vec = field(sim, ω, x_vec, a_vec; basis_order=basis_order)

    # Construct results object
    field_vec = reshape(map(f->SVector{FieldDim,Complex{T}}(f), field_vec), :, 1)
    return FrequencySimulationResult{Dim,FieldDim,P,T}(field_vec, x_vec, RowVector([ω]))

end

function run(sim::FrequencySimulation{Dim,P,T}, ω::Vector{T}, x::Vector{SVector{Dim,T}})::FrequencySimulationResult{Dim,P,T} where {Dim,P,T}
    # Compute for each angular frequency, then join up all the results
    mapreduce(ω->run(sim,ω,x), union, ω)
end

function run(sim::FrequencySimulation{Dim,P,T}, ω::Vector{T}, x::SVector{Dim,T})::FrequencySimulationResult{Dim,P,T} where {Dim,P,T}
    run(sim,ω,[x])
end

function run(sim::FrequencySimulation{Dim,P,T}, ω::T, x::SVector{Dim,T})::FrequencySimulationResult{Dim,P,T} where {Dim,P,T}
    run(sim,[ω],[x])
end

function forcing(source::Source{P,T}, particles::Vector{AbstractParticle{Dim,T}}, t_matrices::Vector{AbstractMatrix}, ω::T, Nh::Integer)::Vector{Complex{T}} where {Dim,P,T}
    mat = [source.coef(n,origin(p),ω) for n in -Nh:Nh, p in particles]
    f = Vector{Complex{T}}(prod(size(mat)))
    H = 2Nh + 1
    for i in eachindex(particles)
        f[((i-1)*H+1):(i*H)] .= t_matrices[i]*mat[:,i]
    end
    return f
end

function basis_coefficients(sim::FrequencySimulation{Dim,P,T}, ω::T; basis_order::Int = 5) where {Dim,P,T}

    # Precompute T-matrices for these particles
    t_matrices = get_t_matrices(sim.medium, sim.particles, ω, basis_order)

    # Compute scattering matrix for all particles
    S = scattering_matrix(sim.medium, sim.particles, t_matrices, ω, basis_order)

    # Get forcing vector for this source
    f = forcing(sim.source, sim.particles, t_matrices, ω, basis_order)

    # Find Hankel coefficients by solving scattering matrix for this forcing
    a = S\f

end

function field(sim::FrequencySimulation{Dim,P,T}, ω::T, x_vec::Vector{SVector{Dim,T}}, a_vec; basis_order::Int=5) where {Dim,P,T}
    Nh = basis_order
    num_particles = length(sim.particles)
    a = OffsetArray(reshape(a_vec,2Nh+1,num_particles),-Nh:Nh,1:num_particles)
    basis = basis_function(sim.medium, ω)
    function sum_basis(x)
        sum(eachindex(sim.particles)) do i
            p = sim.particles[i]
            sum(-Nh:Nh) do m
                a[m,i] * basis(m, x-origin(p))
            end
        end
    end
    map(eachindex(x_vec)) do i
        ind = find(inside(p.shape,x_vec[i]) for p in sim.particles)
        if isempty(ind)
            sim.source.field(x_vec[i],ω) + sum_basis(x_vec[i])
        else
            j = ind[1]
            field(sim.particles[j], sim.medium, ω, x_vec[i], collect(a[:,j]); basis_order=Nh)
        end
    end
end

function field(p::Particle, medium::PhysicalProperties, ω::T, x::AbstractVector{T}, a_vec::AbstractVector; basis_order::Int=5) where T<:Number
    Nh = basis_order

    if inside(p.shape,x)
        inner_basis = basis_function(p, ω)
        b_vec = inner_basis_coefficients(p, medium, ω, a_vec; basis_order=basis_order)
        sum(-Nh:Nh) do m
            inner_basis(m, x-origin(p)) * b_vec[m+Nh+1]
        end
    else
        basis = basis_function(sim.medium, ω)
        sum(-Nh:Nh) do m
            basis(m, x-origin(p)) * a_vec[m+Nh+1]
        end
    end
end

mutable struct TimeSimulation{Dim,P<:PhysicalProperties,T<:AbstractFloat} <: Simulation{Dim,T}
    medium::P
    particles::Vector{AbstractParticle{Dim,T}}
    source::Source{P,T}
    impulse::Function
end
