abstract type Simulation{T,Dim} end

"""
    FrequencySimulation(medium::PhysicalProperties,
                        [particles::AbstractParticles=[],]
                        source::Source)

Build a FrequencySimulation. If particles are not provided, an empty array is used.

After building, you can [`run`](@ref) the simulation to get a [`FrequencySimulationResult`](@ref).
"""
mutable struct FrequencySimulation{T<:AbstractFloat,Dim,P<:PhysicalProperties} <: Simulation{T,Dim}
    "Physical properties of the medium which the whole simulation sits in"
    medium::P
    "Vector of particles, can be of different types"
    particles::AbstractParticles
    "Source"
    source::Source{P,T}
end

# Constructor which infers parametric types from input arguments, note that we  don't need to do much type checking as the struct will error is inconsistent
function FrequencySimulation(medium::P, particles::AbstractParticles{T,Dim}, source::Source{P,T}) where {Dim,T,P<:PhysicalProperties{T,Dim}}
    FrequencySimulation{T,Dim,P}(medium, particles, source)
end

# A simulation with just sources is perfectly reasonable
function FrequencySimulation(medium::P, source::Source{P,T}) where {Dim,T,P<:PhysicalProperties{T,Dim}}
    FrequencySimulation{T,Dim,P}(medium, Vector{AbstractParticle{T,Dim}}(undef,0), source)
end

import Base.show
function show(io::IO, mime::MIME"text/plain", s::FrequencySimulation{T}) where {T}
    # text/plain is repl type
    # FrequencySimulation paramaters can be determined entirely from the medium and shape so we do not need to print them
    println(io, "FrequencySimulation{$T}")
    print(io,   "medium    = ")
    show(io, s.medium)
    println(io)
    println(io, "particles = ")
    for particle in s.particles
        show(io, particle)
    end
    return
end

import Base.run

# Main run function, all other run functions use this
function run(sim::FrequencySimulation{T,Dim,P}, x_vec::Union{Vector{Vector{T}},Vector{SVector{Dim,T}}}, ω::T;
        basis_order::Int = 5) where {Dim,FieldDim,T,P<:PhysicalProperties{T,Dim,FieldDim}}

    x_vec = [SVector{Dim,T}(x...) for x in x_vec]
    # Calculate the Hankel coefficients around each particle, this is where most of the maths happens
    a_vec = basis_coefficients(sim, ω; basis_order=basis_order)

    # Evaluate the total field at the requested x positions
    field_vec = field(sim, ω, x_vec, a_vec)

    # Construct results object
    field_vec = reshape(map(f->SVector{FieldDim,Complex{T}}(f), field_vec), :, 1)
    return FrequencySimulationResult{T,Dim,FieldDim}(field_vec, x_vec, Vector([ω]))

end

function run(sim::FrequencySimulation{T,Dim,P}, x_vec::Union{Vector{Vector{T}},Vector{SVector{Dim,T}}}, ωs::AbstractArray{T}=T[];
        ts::AbstractArray{T} = T[], result_in_time = !isempty(ts),
        basis_order::Int = 5,
        min_basis_order::Int = max(3, Int(round(ωs[1] * basis_order / ωs[end]))),
        basis_order_vec::AbstractVector{Int} = [-1],
        kws...)::(SimulationResult{T,Dim,FieldDim} where FieldDim)  where {Dim,P,T}

    if isempty(ωs) ωs = t_to_ω(ts) end
    x_vec = [SVector{Dim,T}(x...) for x in x_vec]

    # Considering basis_order to be the maximum basis order, then to avoid too much truncation error we use smaller basis orders on the smaller frequencies.
    if basis_order_vec == [-1]
        max_basis_order = max(basis_order,min_basis_order)
        basis_order_vec = Int.(round.(
            LinRange(min_basis_order, max_basis_order, length(ωs))
        ))
        basis_order_vec = basis_order_vec[sortperm(ωs)]
    end

    # ugly bit of code to seperate keywords for simulating frequencies
    ks = []
    freq_kws = Iterators.filter(k -> contains(==,ks,k[1]), kws)
    time_kws = setdiff(kws,freq_kws)

    # if user asks for ω = 0, then we provide
    if first(ωs) == zero(T)
        # Compute for each angular frequency, then join up all the results
        fields = mapreduce(
            i -> run(sim,x_vec,ωs[i]; basis_order = basis_order_vec[i], freq_kws...).field,
        hcat, eachindex(ωs)[2:end])

        # extrapolate field at ω = 0, which should be real when the time signal is real
        zeroresponse = real(ωs[3].*fields[:,1] - ωs[2].*fields[:,2])./(ωs[3]-ωs[2])
        fields = reshape([zeroresponse; fields[:]], length(zeroresponse), size(fields,2)+1)
    else
        fields = mapreduce(
            i -> run(sim,x_vec,ωs[i]; basis_order = basis_order_vec[i], freq_kws...).field,
        hcat, eachindex(ωs))
    end

    if !result_in_time
        FrequencySimulationResult(fields,x_vec,Vector(ωs))
    else
        if isempty(ts) ts = ω_to_t(ωs) end

        # better to use the defaults of TimeSimulationResult's Constructor.
        frequency_to_time(FrequencySimulationResult(fields,x_vec,Vector(ωs)); t_vec = reshape(ts,length(ts)), time_kws...)
    end
end

function run(sim::FrequencySimulation{T,Dim,P}, x::AbstractVector{T}, ωs::AbstractVector{T}=T[];
        kws...)::(SimulationResult{T,Dim,FieldDim} where FieldDim) where {Dim,P,T}
    run(sim,[x],ωs; kws...)
end

function run(sim::FrequencySimulation{T,Dim,P}, x::AbstractVector{T}, ω::T;
        kws...)::(SimulationResult{T,Dim,FieldDim} where FieldDim) where {Dim,P,T}
    run(sim,[x],[ω]; kws...)
end

# Add docstring to run functions
"""
    run(sim::FrequencySimulation, x, ω; basis_order=5)

Run the simulation `sim` for the position `x` and angular frequency `ω`.

Position can be an SVector or Vector{SVector} and frequency can be a float or
vector of floats.
"""
function run(s::FrequencySimulation) throw(MethodError(run, (s,))) end

"""
    run(sim::FrequencySimulation, rectangle;
        res=20, xres=res, yres=res, basis_order=5)

Run the simulation `sim` for a grid of positions in rectangle and for angular frequency `ω`.

Frequency can be a float or vector of floats. The resolution of the grid points is defined
by xres and yres.
"""
function run(sim::FrequencySimulation, rect::Rectangle, ω_vec::AbstractVector;
             res=20, xres=res, yres=res, kws...)

    #Size of the step in x and y direction
    step_size = [rect.width / xres, rect.height / yres]
    x_vec = [SVector(bottomleft(rect) + step_size.*[i,j]) for i=0:xres, j=0:yres]

    return run(sim, x_vec[:], ω_vec; kws...)
end

"""
    forcing(source::Source, particles::AbstractParticles, ω::AbstractFloat, Nh::Integer)::Vector{Complex}

Create forcing vector from source, forms the right hand side of matrix equation to find [`basis_coefficients`](@ref).
"""
function forcing(source::Source{P,T}, particles::AbstractParticles, ω::T, Nh::Integer)::Vector{Complex{T}} where {P,T}
    mat = [source.coef(n,origin(p),ω) for n in -Nh:Nh, p in particles]
    f = Vector{Complex{T}}(undef,prod(size(mat)))
    H = 2Nh + 1
    for i in eachindex(particles)
        f[((i-1)*H+1):(i*H)] .= mat[:,i]
    end
    return f
end

"""
    basis_coefficients(sim::FrequencySimulation, ω::AbstractFloat; basis_order::Int=5)::Matrix{Complex}

Return coefficients for bases around each particle for a given simulation and angular frequency (ω).
"""
function basis_coefficients(sim::FrequencySimulation{T,Dim,P}, ω::T; basis_order::Int = 5) where {Dim,P,T}

    # Precompute T-matrices for these particles
    t_matrices = get_t_matrices(sim.medium, sim.particles, ω, basis_order)

    # Compute scattering matrix for all particles
    S = scattering_matrix(sim.medium, sim.particles, t_matrices, ω, basis_order)

    # Get forcing vector for this source
    f = forcing(sim.source, sim.particles, ω, basis_order)

    # Find Hankel coefficients by solving scattering matrix for this forcing
    a = S\f

    # reshape and multiply by t-matrix to get the scattering coefficients
    a = reshape(a,2basis_order+1,length(sim.particles))
    for i in axes(a,2)
        a[:,i] = t_matrices[i] * a[:,i]
    end
    a
end

function field(sim::FrequencySimulation{T,Dim,P}, ω::T, x_vec::Vector{SVector{Dim,T}}, a_vec) where {Dim,P,T}

    Nh = Int((size(a_vec,1) - one(T)) / T(2.0)) # basis_order
    num_particles = length(sim.particles)
    a = OffsetArray(a_vec,-Nh:Nh,1:num_particles)
    basis = basis_function(sim.medium, ω)

    function sum_basis(x)
        sum(eachindex(sim.particles)) do i
            p = sim.particles[i]
            sum(-Nh:Nh) do m
                a[m,i] * basis(m, x-origin(p))
            end
        end
    end
    map(x_vec) do x
        j = findfirst(p -> x∈p, sim.particles)
        if typeof(j) === Nothing
            sim.source.field(x,ω) + (isempty(sim.particles) ? zero(Complex{T}) : sum_basis(x))
        else
            p = sim.particles[j]
            internal_field(x, p, sim, ω, collect(a[:,j]))
        end
    end
end
