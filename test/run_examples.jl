@testset "run examples and generate README.md" begin
    examplestoconvert = ["random_particles", "particles_in_circle"]

    for str in examplestoconvert
        include(string("../example/",str,"/",str,".jl"))
        @test true
    end

    # auto generate their README.md to guarantee the code in these READMEs works.

    include("../example/generate-READMEs.jl")

    @test true
end
