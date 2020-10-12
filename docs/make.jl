using Documenter, MultipleScattering

makedocs(
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ), # changes urls if built locally
    sitename = "MultipleScattering.jl",
    authors = "Artur L. Gower and Jonathan Deakin",
    source = "src",
    # modules = [MultipleScattering],
    pages = [
        "Home" =>"index.md",
        "Manual" => [
            "manual/intro.md",
            "manual/source.md",
            "manual/shapes.md",
            "manual/particles.md",
            "manual/time_response.md",
            "manual/plot.md",
            "manual/new_types.md"
        ],
        "Library" => [
            "library/base.md",
            "library/acoustics.md",
            "library/random.md"
        ],
        "Theory" => "maths/README.md",
        "Examples" =>
        [
            "example/README.md",
            "example/box_size/README.md",
            "example/hankel_convergence/README.md",
            "example/moments/README.md",
            "example/near_surface_backscattering/README.md",
            "example/random_particles/README.md",
            "example/time_response_single_particle/README.md",
            "example/two_particles/README.md"
        ]
    ]
)

deploydocs(
    branch = "gh-pages",
    target = "build",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#.#", devurl => devurl],
    repo = "github.com/JuliaWaveScattering/MultipleScattering.jl.git"
)
