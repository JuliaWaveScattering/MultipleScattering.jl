using MultipleScattering, Test

@testset "run examples and generate README.md" begin
    examplestoconvert = ["random_particles", "particles_in_circle"]
    folder = "../docs/src/example/"

    for str in examplestoconvert
        include(string(folder,str,"/",str,".jl"))
        @test true
    end

    # auto generate their README.md to guarantee the code in these READMEs works.

    include(folder*"generate-READMEs.jl")

    gernerate_README(examplestoconvert, folder)

    @test true
end
