using Documenter, MultipleScattering

makedocs(
    format=Documenter.HTML(),
    sitename="MultipleScattering.jl",
    authors = "Artur Gower and Jonathan Deakin",
    source= "src",
    modules=[MultipleScattering],
    pages=[
        "Home" =>"index.md",
        "Manual" => [
            "base.md",
            "acoustics.md",
            "random.md"
        ],
        "Examples" => [
            "example/box_size/README.md",
            "example/hankel_convergence/README.md",
            "example/intro/README.md",
            "example/lens/README.md",
            "example/moments/README.md",
            "example/near_surface_backscattering/README.md",
            "example/random_particles/README.md",
            "example/random_particles_in_circle/README.md",
            "example/shapes/README.md",
            "example/time_response_single_particle/README.md",
            "example/two_particles/README.md",
        ]
    ]
)

deploydocs(
    branch = "gh-pages",
    latest = "master",
    julia = "1.1",
    osname = "linux",
    target = "build",
    deps = nothing,
    make = nothing,
    repo = "github.com/jondea/MultipleScattering.jl.git"
)
