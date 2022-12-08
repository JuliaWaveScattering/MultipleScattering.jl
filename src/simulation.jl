abstract type Simulation{Dim} end

spatial_dimension(::Simulation{Dim}) where {Dim} = Dim

"""
    FrequencySimulation([particles::AbstractParticles=[],]
                        source::AbstractSource)

Build a FrequencySimulation. If particles are not provided, an empty array is used.

After building, you can [`run`](@ref) the simulation to get a [`FrequencySimulationResult`](@ref).
"""
mutable struct FrequencySimulation{Dim,P<:PhysicalMedium} <: Simulation{Dim}
    "Vector of particles, can be of different types."
    particles::AbstractParticles
    "RegularSource wave, where source.medium is the background medium of the simulation."
    source::AbstractSource{P}
end

field_dimension(::FrequencySimulation{Dim,P}) where {Dim,P} = field_dimension(P)

# Constructor which infers parametric types from input arguments, note that we  don't need to do much type checking as the struct will error is inconsistent
function FrequencySimulation(particles::AbstractParticles{Dim}, source::AbstractSource{P}) where {Dim,P<:PhysicalMedium{Dim}}
    FrequencySimulation{Dim,P}(particles, source)
end

# A simulation with just sources is perfectly reasonable
function FrequencySimulation(source::AbstractSource{P}) where {Dim,P<:PhysicalMedium{Dim}}
    FrequencySimulation{Dim,P}(AbstractParticle{Dim}[], source)
end

import Base.show
function show(io::IO, mime::MIME"text/plain", s::FrequencySimulation{T}) where {T}
    # text/plain is repl type
    # FrequencySimulation paramaters can be determined entirely from the medium and shape so we do not need to print them
    println(io, "FrequencySimulation{$T}")
    print(io,   "medium    = ")
    show(io, s.source.medium)
    println(io)
    println(io, "particles = ")
    for particle in s.particles
        show(io, particle)
    end
    return
end

import Base.run

run(particles::AbstractParticles, source::AbstractSource, x::Union{Shape,AbstractVector}, ω; kws...) = run(FrequencySimulation(particles,source), x, ω; kws...)

run(source::AbstractSource, x::Union{Shape,Vector}, ω; kws...) = run(FrequencySimulation(source), x, ω; kws...)

function run(sim::FrequencySimulation, x::AbstractVector{<:Number}, ωs::AbstractVector{<:Number}=eltype(x)[];
        kws...)::SimulationResult
    run(sim,[x],ωs; kws...)
end

function run(sim::FrequencySimulation, x::AbstractVector{<:Number}, ω::Number;
        kws...)::SimulationResult
    run(sim,[x],[ω]; kws...)
end

# Main run function, all other run functions use this
function run(sim::FrequencySimulation, x_vec::Vector{V}, ω::Number;
        basis_order::Integer = 5,
        only_scattered_waves::Bool = false, kws...) where V<:AbstractVector

    Dim = spatial_dimension(sim)
    FieldDim = field_dimension(sim)

    x_vec = [SVector{Dim}(x) for x in x_vec]

    # return just the source if there are no particles
    if isempty(sim.particles)
        f = field(sim.source)
        field_vec = f.(x_vec,ω)
    else
        # Calculate the outgoing basis coefficients around each particle, this is where most of the maths happens
        a_vec = basis_coefficients(sim, ω; basis_order=basis_order)

        if only_scattered_waves # remove source wave
            sim = FrequencySimulation(sim.particles,
                constant_source(sim.source.medium, 0)
            )
        end

        # Evaluate the total field at the requested x positions
        field_vec = field(sim, ω, x_vec, a_vec)
    end

    # Construct results object
    T = real(eltype(field_vec))
    field_vec = reshape(map(f->SVector{FieldDim,Complex{T}}(f), field_vec), :, 1)
    return FrequencySimulationResult{T,Dim,FieldDim}(field_vec, x_vec, Vector([ω]))

end

function run(sim::FrequencySimulation, x_vec::Vector{V}, ωs::AbstractArray=[];
        ts::AbstractArray = [], result_in_time = !isempty(ts),
        basis_order::Int = 5,
        min_basis_order::Int = basis_order,
        basis_order_vec::AbstractVector{Int} = [-1],
        kws...)::SimulationResult where V<:AbstractVector

    Dim = spatial_dimension(sim)

    if isempty(ωs) ωs = t_to_ω(ts) end
    x_vec = [SVector{Dim}(x) for x in x_vec]

    # Considering basis_order to be the maximum basis order, then to avoid too much truncation error we use smaller basis orders on the smaller frequencies.
    if basis_order_vec == [-1]
        max_basis_order = max(basis_order,min_basis_order)
        basis_order_vec = Int.(round.(
            LinRange(min_basis_order, max_basis_order, length(ωs))
        ))
        basis_order_vec = basis_order_vec[sortperm(ωs)]
    end

    # if user asks for ω = 0, then we provide
    if first(ωs) == 0
        # Compute for each angular frequency, then join up all the results
        fields = mapreduce(
            i -> run(sim,x_vec,ωs[i]; basis_order = basis_order_vec[i], kws...).field,
        hcat, eachindex(ωs)[2:end])

        # extrapolate field at ω = 0, which should be real when the time signal is real
        zeroresponse = real(ωs[3].*fields[:,1] - ωs[2].*fields[:,2]) ./ (ωs[3]-ωs[2])
        fields = reshape([zeroresponse; fields[:]], length(zeroresponse), size(fields,2)+1)
    else
        fields = mapreduce(
            i -> run(sim,x_vec,ωs[i]; basis_order = basis_order_vec[i], kws...).field,
        hcat, eachindex(ωs))
    end

    if !result_in_time
        FrequencySimulationResult(fields,x_vec,Vector(ωs))
    else
        @error "Using `result_in_time=true` in `run` is depricated and will be removed in later versions. Please instead use `frequency_to_time` on the results of `run`."
        # if isempty(ts) ts = ω_to_t(ωs) end
        #
        # # better to use the defaults of TimeSimulationResult's Constructor.
        # frequency_to_time(FrequencySimulationResult(fields,x_vec,Vector(ωs)); t_vec = reshape(ts,length(ts)), time_kws...)
    end
end

# Add docstring to run functions
"""
    run(sim::FrequencySimulation, x, ω; basis_order=5)

Run the simulation `sim` for the position `x` and angular frequency `ω`.

Position can be an SVector or Vector{SVector} and frequency can be a float or
vector of floats.
"""
function run(s::FrequencySimulation) throw(MethodError(run, (s,))) end

run(sim::FrequencySimulation, region::Shape, ω::AbstractFloat; kws...) = run(sim, region, [ω]; kws...)

"""
    run(sim::FrequencySimulation, region::Shape;
        res=20, xres=res, yres=res, basis_order=5)

Run the simulation `sim` for a grid of positions in region and for angular frequency `ω`.

Frequency can be a float or vector of floats. The resolution of the grid points is defined
by xres and yres.
"""
function run(sim::FrequencySimulation, region::Shape, ω_vec::AbstractVector;
            kws...)

    x_vec, inds = points_in_shape(region; kws...)
    # x_vec is a square grid of points and x_vec[inds] are the points in the region.

    result = run(sim, x_vec[inds], ω_vec; kws...)

    field_mat = zeros(typeof(result.field[1]), length(x_vec), length(ω_vec))
    field_mat[inds,:] = result.field

    return FrequencySimulationResult(field_mat, x_vec, ω_vec)
end

"""
    basis_coefficients(sim::FrequencySimulation, ω::AbstractFloat; basis_order::Int=5)::Matrix{Complex}

Return coefficients for bases around each particle for a given simulation and angular frequency (ω).
"""
function basis_coefficients(sim::FrequencySimulation{Dim,P}, ω::Number; basis_order::Int = 5) where {Dim,P}

    # Precompute T-matrices for these particles
    t_matrices = get_t_matrices(sim.source.medium, sim.particles, ω, basis_order)

    # Compute scattering matrix for all particles
    S = scattering_matrix(sim.source.medium, sim.particles, t_matrices, ω, basis_order)

    # Get forcing vector from source, forms the right hand side of matrix equation to find basis_coefficients
    source_coefficient = regular_spherical_coefficients(sim.source)
    forcing = reduce(vcat, [source_coefficient(basis_order,origin(p),ω) for p in sim.particles])

    # Find scattering coefficients by solving this forcing
    a = (S + I) \ forcing

    # reshape and multiply by t-matrix to get the scattering coefficients
    a = reshape(a,basisorder_to_basislength(P,basis_order),length(sim.particles))
    for i in axes(a,2)
        a[:,i] = t_matrices[i] * a[:,i]
    end
    return a
end

function field(sim::FrequencySimulation{Dim,P}, ω::Number, x_vec::Vector{V}, a_vec) where {Dim,P,T<:Number,V <: AbstractArray{T}}

    N = basislength_to_basisorder(P,size(a_vec,1))
    num_particles = length(sim.particles)
    basis = outgoing_basis_function(sim.source.medium, ω)

    function sum_basis(x)
        sum(eachindex(sim.particles)) do i
            p = sim.particles[i]
            sum(a_vec[:,i] .* basis(N, x-origin(p)))
        end
    end
    map(x_vec) do x
        j = findfirst(p -> x ∈ p, sim.particles)
        if typeof(j) === Nothing
            field(sim.source)(x,ω) + sum_basis(x)
        else
            p = sim.particles[j]
            internal_field(x, p, sim.source, ω, collect(a_vec[:,j]))
        end
    end
end
