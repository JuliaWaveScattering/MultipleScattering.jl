# Simulations and results

```@meta
DocTestSetup = quote
    using MultipleScattering
end
```

This page explains what happens when you [`run`](@ref) a [`FrequencySimulation`](@ref), which keyword arguments control the accuracy, and how to use the results. For a first example see the [Introduction](@ref).

## Setting up a simulation

A [`FrequencySimulation`](@ref) needs a vector of particles and a source, where the medium of the source is the background (host) medium:

```jldoctest sim
host_medium = Acoustic(2; ρ = 1.0, c = 1.0)
particle_medium = Acoustic(2; ρ = 0.5, c = 0.6)

particles = [
    Particle(particle_medium, Sphere([0.0, 0.0], 1.0)),
    Particle(particle_medium, Sphere([3.0, 1.0], 0.5))
]

source = plane_source(host_medium; direction = [1.0, 0.0])
sim = FrequencySimulation(particles, source)

# output

FrequencySimulation{2}
medium    = Acoustic(1.0, 1.0 + 0.0im, 2)
particles = 
Particle(Acoustic(0.5, 0.6 + 0.0im, 2), Sphere{Float64, 2}([0.0, 0.0], 1.0))
Particle(Acoustic(0.5, 0.6 + 0.0im, 2), Sphere{Float64, 2}([3.0, 1.0], 0.5))
```

A simulation with no particles, `FrequencySimulation(source)`, is also allowed and just returns the source field. Particles are not allowed to overlap; [`overlapping_pairs`](@ref) returns the pairs of particles which do.

## [What `run` calculates](@id run_maths)

Let $\mathbf x_j$ be the origin of particle $j$, $k = \omega / c$ the wavenumber of the host medium, and let $\mathrm v_n$ and $\mathrm u_n$ be the regular and outgoing basis functions, given by [`regular_basis_function`](@ref) and [`outgoing_basis_function`](@ref). In 2D

```math
\mathrm v_n(\mathbf x) = \mathrm J_n(k r) \mathrm e^{\mathrm i n \theta}, \qquad \mathrm u_n(\mathbf x) = \mathrm H^{(1)}_n(k r) \mathrm e^{\mathrm i n \theta},
```

where $(r, \theta)$ are polar coordinates, and in 3D $\mathrm j_\ell(k r) \mathrm Y_\ell^m(\hat{\mathbf x})$ and $\mathrm h^{(1)}_\ell(k r) \mathrm Y_\ell^m(\hat{\mathbf x})$, with spherical Bessel and Hankel functions and spherical harmonics.

Outside of the particles, the total field is the source field plus the waves scattered by each particle:

```math
u(\mathbf x) = u_{\text{in}}(\mathbf x) + \sum_{j} \sum_{n} f^j_n \, \mathrm u_n(\mathbf x - \mathbf x_j).
```

To find the scattering coefficients $f^j_n$, `run` does the following:

1. Calculates the T-matrix $T^j$ of each particle with [`t_matrix`](@ref). This matrix gives the waves scattered by a particle in terms of the waves incident on it. Congruent particles share the same T-matrix, see [`get_t_matrices`](@ref).
2. Expands the source around each particle, $u_{\text{in}}(\mathbf x) = \sum_n a^j_n \mathrm v_n(\mathbf x - \mathbf x_j)$, where the vector of coefficients $a^j$ is given by `source.coefficients(basis_order, x_j, ω)`.
3. Solves the multiple scattering system: the wave incident on particle $j$ is the source plus the waves scattered by all other particles,
   ```math
   f^j = T^j g^j, \qquad g^j = a^j + \sum_{l \not = j} U(\mathbf x_j - \mathbf x_l)^{\mathsf T} f^l,
   ```
   where $U$ is given by [`outgoing_translation_matrix`](@ref). This is the linear system built with [`scattering_matrix`](@ref) and solved in [`basis_coefficients`](@ref).
4. Evaluates the field at the requested points. Inside a particle it uses [`internal_field`](@ref) instead.

The coefficients $f^j_n$ are returned by `basis_coefficients`, with one column for each particle:

```jldoctest sim
julia> ω = 2.0;

julia> f = basis_coefficients(sim, ω; basis_order = 4);

julia> size(f) # 2*4 + 1 = 9 coefficients in 2D for each particle
(9, 2)
```

For the details of this method see the notes in [Theory](@ref).

## Choosing the basis order

The sums over $n$ above are truncated at the basis order $N$, using $|n| \leq N$ in 2D and $\ell \leq N$ in 3D. This is the most important parameter for the accuracy and the speed of a simulation. The default is `basis_order = 5`.

A larger basis order is needed for larger particles or higher frequencies. A good first guess is given by [`estimate_outgoing_basisorder`](@ref), which checks the convergence of the T-matrix of a particle:

```jldoctest sim
julia> estimate_outgoing_basisorder(host_medium, particles[1], ω)
4
```

When particles are close together the waves between them vary more rapidly, so it is good practice to check that the results do not change when increasing the basis order:

```jldoctest sim
julia> x = [5.0, 0.0];

julia> us = [field(run(sim, x, ω; basis_order = N))[1] for N in [2, 4, 8, 12]];

julia> round.(abs.(us[1:3] .- us[4]), sigdigits = 2) # difference from using basis_order = 12
3-element Vector{Float64}:
 0.0038
 0.00012
 4.3e-11
```

## Running for many positions and frequencies

The positions can be one vector or a vector of vectors, and the angular frequencies can be one number or a vector:

```jldoctest sim
julia> result = run(sim, [[5.0, 0.0], [0.0, 3.0]], [1.0, 2.0]);

julia> size(result) # (number of positions, number of frequencies)
(2, 2)

julia> round.(field(result), digits = 4)
2×2 Matrix{ComplexF64}:
 0.6596-0.4263im  0.2426-0.681im
 0.9771-0.1341im  1.1626-0.5013im
```

You can also skip creating a `FrequencySimulation` and use `run(particles, source, x, ωs)`.

When running for many frequencies, the keyword `basis_order` is used for the largest frequency. To save time, smaller basis orders can be used for smaller frequencies with `min_basis_order`. Then the basis order increases linearly from `min_basis_order`, for the smallest frequency, up to `basis_order`, for the largest. Alternatively, give one basis order for each frequency with `basis_order_vec`.

If the first frequency is zero, then the field for `ω = 0` is extrapolated from the next two frequencies, as the T-matrices are often not defined at `ω = 0`. This is useful when converting to the [Time response](@ref).

## Running over a region

To calculate the field over a region, pass a [`Shape`](@ref) instead of positions:

```jldoctest sim
julia> region = Box([[-2.0, -2.0], [5.0, 3.0]]);

julia> result = run(sim, region, ω; res = 10);

julia> size(result) # an 11 x 11 grid of points
(121, 1)
```

The grid of points covers the [`bounding_box`](@ref) of the region, and the field is set to zero for points which are not inside the region. Use `xres` and `yres` to choose the resolution along each axis, and `exclude_region` to avoid evaluating the field in some shape, for example close to a point source. In 3D the field is calculated on the $x$-$z$ plane at the height `y`, which defaults to `y = 0`, using `xres` and `zres`. See [`points_in_shape`](@ref) for details.

## Only the scattered waves

Use `only_scattered_waves = true` to remove the source from the field outside of the particles:

```jldoctest sim
julia> x_out = [5.0, 0.0];

julia> u_total = field(run(sim, x_out, ω; basis_order = 8))[1];

julia> u_scat = field(run(sim, x_out, ω; basis_order = 8, only_scattered_waves = true))[1];

julia> u_total ≈ u_scat + field(source, x_out, ω)
true
```

Inside the particles the total field is always returned.

## Results

`run` returns a [`FrequencySimulationResult`](@ref), with the fields:

- `result.x`: the positions, a vector of `SVector`;
- `result.ω`: the angular frequencies;
- `result.field`: a matrix of `SVector` with one row for each position and one column for each frequency.

For scalar waves, such as acoustics, [`field`](@ref) returns a matrix of complex numbers instead, and `field(result, i, j)` returns the field at position `i` and frequency `j`. Results can be added together, or multiplied by numbers, for example `2.0 * result1 + result2`.

To convert the results to the time domain, use [`frequency_to_time`](@ref), see [Time response](@ref). To plot the results see [Plotting](@ref).
