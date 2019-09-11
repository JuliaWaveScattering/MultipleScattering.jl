# Source wave

```@meta
DocTestSetup = quote
    using MultipleScattering
end
```
[`Source`](@ref) is a `struct` a For acoustics, any wave field $u(x,y)$ that satisfies $\nabla^2 u(x,y) + k^2 u(x,y) = 0$, with $k = \omega/c$, can be a source wave, also called an incident wave. See [Source](@ref) for a list of types and functions.

## 2D Acoustics

Two common source waves are shown below.

For a plane-wave of the form $u(x,y) = A \mathrm e^{\mathrm i k \mathbf n \cdot (\mathbf x - \mathbf x_0)}$, where $A$ is the amplitude, $\mathbf n = (n_1,n_2)$ is the direction of propagation, and $\mathbf x_0 = (x_0,y_0)$ is the initially position of the source, we can use
```jldoctest intro
julia> medium = Acoustic(1.0, 1.0, 2);

julia> A = 1.0;

julia> n = [1.0,1.0];

julia> x0 = [1.0,0.0];

julia> source = plane_source(medium; amplitude = A, direction = n, position = x0);
```

We can plot this source wave field for one frequency ω by using
```julia
julia> ω = 1.0;

julia> plot(sim, ω; bounds = Rectangle([-1.0,-1.0],[1.0,1.0]))
```

plot(sim::FrequencySimulation{T}, ω::T;

source_field(x,ω) = amp(ω)*exp(im*ω/medium.c*dot(x-position, direction))


source_position = SVector(0.0,1.0)
amplitude = 1.0
s1 = point_source(a2, source_position, amplitude)
s2 = point_source(a2, 2.0*source_position, amplitude)

source = plane_source(host_medium; direction = [1.0,0.0])
