using Documenter, MultipleScattering

makedocs(
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ), # changes urls if built locally
    sitename="MultipleScattering.jl",
    authors = "Artur L. Gower and Jonathan Deakin",
    source= "src",
    modules=[MultipleScattering],
    pages=[
        "Home" =>"index.md",
        "Manual" => [
            "manual/intro.md",
            "manual/shapes.md",
            "manual/time_response.md",
            "manual/plot.md"
        ],
        "Library" => [
            "library/base.md",
            "library/acoustics.md",
            "library/random.md"
        ],
        "Examples" =>
        [
            "example/README.md",
            "example/box_size/README.md",
            "example/hankel_convergence/README.md",
            "example/intro/README.md",
            "example/lens/README.md",
            "example/moments/README.md",
            "example/near_surface_backscattering/README.md",
            "example/random_particles/README.md",
            "example/random_particles_in_circle/README.md",
            "example/time_response_single_particle/README.md",
            "example/two_particles/README.md",
        ]
    ]
)

deploydocs(
    branch = "gh-pages",
    latest = "master",
    target = "build",
    repo = "github.com/jondea/MultipleScattering.jl.git"
)
