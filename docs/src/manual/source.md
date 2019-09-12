# Source wave

```@meta
DocTestSetup = quote
    using MultipleScattering
end
```
[`Source`](@ref) is a `struct` a For acoustics, any wave field $u(x,y)$ that satisfies $\nabla^2 u(x,y) + k^2 u(x,y) = 0$, with $k = \omega/c$, can be a source wave, also called an incident wave. See [`Source`](@ref) for a list of types and functions.

## 2D Acoustics

Two common source waves are shown below.

For a plane-wave of the form $u(x,y) = A \mathrm e^{\mathrm i k \mathbf n \cdot (\mathbf x - \mathbf x_0)}$, where $A$ is the amplitude, $\mathbf n = (n_1,n_2)$ is the direction of propagation, and $\mathbf x_0 = (x_0,y_0)$ is the initially position of the source, we can use
```jldoctest intro
julia> medium = Acoustic(1.0, 1.0, 2);

julia> A = 1.0;

julia> n = [1.0,1.0];

julia> x0 = [1.0,0.0];

julia> plane_wave = plane_source(medium; amplitude = A, direction = n, position = x0);
```

We can plot this source wave one frequency ω by using
```julia
julia> ω = 1.0;

julia> domain = Rectangle([-1.0,-1.0],[1.0,1.0]);

julia> plot(plane_wave, ω; bounds = domain)
```
![Plot plane wave](../media/plane-wave.png)

Another useful source is the point source $u(x,y) = \frac{\mathrm i A}{4} \mathrm H_0^{(1)}(k \|(x-x_0,y-y_0)\|)$ where $A$ is the amplitude,  $\mathbf x_0 = (x_0,y_0)$ is the origin of the point source, and $\mathrm H_0^{(1)}$ is the Hankel function of the first kind.

```jldoctest intro
julia> x0 = [0.0,-1.2];

julia> domain = Rectangle([-1.0,-1.0],[1.0,1.0]);

julia> point_wave = point_source(medium, x0, A);
```
```julia
julia> plot(point_wave, ω; bounds = domain)
```
![Plot point wave](../media/point-wave.png)

!!! note
    Because the point source has a singularity at $x_0$ it is best to avoid plotting, and evaluating the field, close to $x_0$.

## Creating new sources

The easiest way to create new sources is to just sum together predefined sources:

```julia
julia> source = (3.0 + 1.0im) * point_wave + plane_wave;

julia> plot(source, ω; bounds = domain)
```
![Plot point wave](../media/combined-source.png)

For example, to create a finite emitter/transducer source we can use:
```julia
julia> ys = LinRange(-0.7, 0.7, 30);

julia> source = sum(ys) do y point_source(medium, [-1.1, y]) end;

julia> plot(source, 4.0; bounds = domain, field_apply = abs, res = 40)
```
![Plot point wave](../media/transducer-source.png)

where `field_apply` is applied to the wave field at every point, the default is `field_apply = real`, and `res` is the resolution along both the $x$ and $y$ axis.

To define a new source you will need to understand the internals below.

## Source internals

The `struc` [`Source`](@ref) has three fields, shown with examples below:
```jldoctest intro
julia> plane_wave = plane_source(Acoustic(1.0, 1.0, 2); direction = [1.0,0.0]);

julia> plane_wave.medium # the physical medium
Acoustic(1.0, 1.0 + 0.0im, 2)

julia> x = [1.0,1.0]; ω = 1.0;

julia> plane_wave.field(x,ω) # the value of the field
0.5403023058681398 + 0.8414709848078965im
```
However, to calculate the scattering from a particle due to a source, we need to represent the source in a radial coordinate system. That is,  



Source{P<:PhysicalProperties,T<:AbstractFloat}

 `Source.medium` is the physical medium, `Source.field` is a function of`(x,ω)` which gives the source field, `Soruce.coef` is a function of `(x,ω)` which gives the source field expanded in terms of
