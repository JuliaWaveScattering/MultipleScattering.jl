@testset "acoustic sources" begin
    a2 = Acoustic(0.1,0.1 + 0.0im,2)
    a2_host = Acoustic(1.3,1.5 + 0.0im,2)

    # Create two point sources
    source_position = [0.0,1.0]
    amplitude = 1.0
    s1 = point_source(a2, source_position, amplitude)
    s2 = point_source(a2, 2.0*source_position, amplitude)

    # Create new souce as a linear combination of two other sources
    s3 = 2*s1 + s2;

    # Check that the field is indeed a linear conbination
    x = [1.0, 1.0]
    @test field(s3,x,1.0) == 2*field(s1,x,1.0) + field(s2,x,1.0)

    # Test the bessel expansions of the source
    ω = 0.8
    basis_order = 7
    centre =  [1.0,0.0]

    s3_expand = source_expand(s3, centre; basis_order = 7)

    xs = [centre + 0.1 .* [cos(τ),sin(τ)] for τ = 0.0:0.3:1.5]
    @test norm([field(s3,x,ω) - s3_expand(x,ω) for x in xs]) < 1e-7*norm([field(s3,x,ω) for x in xs])

    source = plane_source(a2_host; position = SVector(-10.0,0.0), direction = SVector(1.0,0.0), amplitude = 1.0)
    s_expand = source_expand(source, centre; basis_order = 4)
    @test norm([field(source,x,ω) - s_expand(x,ω) for x in xs]) < 2e-9*norm([field(source,x,ω) for x in xs])


    # Test equivalence between plane-Sources
    # a plane-wave source with explicit wave direction

    direction = SVector(-10.0,1.0)
    position = SVector(rand(2)...)
    amplitude = rand() + rand() * im

    a3_host = Acoustic(3, ρ = 1.3, c = 1.5 + 0.0im)
    @test_throws(DimensionMismatch,PlaneSource(a3_host; direction = direction))
    @test_throws(DimensionMismatch,PlaneSource(a3_host; amplitude = [1.0,2.0]))

    psource = PlaneSource(a2_host; direction = direction, position =  position, amplitude = amplitude)

    source = plane_source(a2_host; direction = direction, position =  position, amplitude = amplitude)

    pf = field(psource)
    sf = field(source)

    xs = [rand(2) for i = 1:10]
    ωs = [rand() for i = 1:10]

    errs = map(xs,ωs) do x,ω
        abs(pf(x,ω) - sf(x,ω))
    end

    @test maximum(errs) ≈ 0.0
end
