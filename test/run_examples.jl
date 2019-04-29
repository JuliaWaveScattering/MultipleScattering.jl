examplestoconvert = ["random_particles", "particles_in_circle"]

for str in examplestoconvert
    include(string("../example/",str,"/",str,".jl"))
    @test true
end
