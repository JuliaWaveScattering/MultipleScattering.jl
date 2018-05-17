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

function field(sim, ω, x_vec, a_vec; basis_order::Int=5)
    num_particles = length(sim.particles)
    a = OffsetArray(reshape(a_vec,2basis_order+1,num_particles),-basis_order:basis_order,1:num_particles)
    basis = basis_function(sim.medium, ω)
    function sum_hankel_coefficients(x)
        sum(eachindex(sim.particles)) do i
            p = sim.particles[i]
            sum(-basis_order:basis_order) do m
                a[m,i] * basis(m, x-origin(p))
            end
        end
    end
    map(sum_hankel_coefficients, x_vec)
end

mutable struct TimeSimulation{Dim,P<:PhysicalProperties,T<:AbstractFloat} <: Simulation{Dim,T}
    medium::P
    particles::Vector{AbstractParticle{Dim,T}}
    source::Source{P,T}
    impulse::Function
end
