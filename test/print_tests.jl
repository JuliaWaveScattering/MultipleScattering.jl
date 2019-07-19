@testset "Print/show" begin
    # Test whether we can reconstruct the objects from their show values

    "Simulates a user printing an object in the REPL, copying and pasting the result and pressing enter"
    function show_and_eval(object)
        io = IOBuffer()
        show(io, object)
        io |> take! |> String |> Meta.parse |> eval
    end

    acoustic = Acoustic(1.89, 1.5, 2)
    circle = Circle(2.0)
    particle = Particle(acoustic, circle)
    simulation = FrequencySimulation(acoustic, [particle], plane_source(acoustic))

    @test acoustic == show_and_eval(acoustic)
    @test circle == show_and_eval(circle)
    @test particle == show_and_eval(particle)

    # Currently can't do this with FrequencySimulations because of Source type
    @test show(devnull, simulation) == nothing

end
