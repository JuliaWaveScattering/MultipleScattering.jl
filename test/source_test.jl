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

    # Test regular basis expansion of the source
    ω = 0.8
    basis_order = 7
    centre =  [1.0,0.0]

    s3_expand = source_expand(s3, centre; basis_order = 7)

    xs = [centre + 0.1 .* [cos(τ),sin(τ)] for τ = 0.0:0.3:1.5]
    @test norm([field(s3,x,ω) - s3_expand(x,ω) for x in xs]) < 1e-7*norm([field(s3,x,ω) for x in xs])

    source = plane_source(a2_host; position = SVector(-10.0,0.0), direction = SVector(1.0,0.0), amplitude = 1.0)
    s_expand = source_expand(source, centre; basis_order = 4)
    @test norm([field(source,x,ω) - s_expand(x,ω) for x in xs]) < 2e-9*norm([field(source,x,ω) for x in xs])

    # Test DimensionMismatch between 3D and 2D

    direction = SVector(-10.0,1.0)
    position = SVector(rand(2)...)
    amplitude = rand() + rand() * im

    a3_host = Acoustic(3, ρ = 1.3, c = 1.5 + 0.0im)
    @test_throws(DimensionMismatch,PlaneSource(a3_host; direction = direction))
    @test_throws(DimensionMismatch,PlaneSource(a3_host; amplitude = [1.0,2.0]))

    # Test equivalence between plane-Sources
    psource = PlaneSource(a2_host; direction = direction, position =  position, amplitude = amplitude)

    source = plane_source(a2_host; direction = direction, position =  position, amplitude = amplitude)

    pf = field(psource)
    sf = field(source)

    xs = [rand(2) for i = 1:10]
    ωs = [rand() for i = 1:10]

    errs = map(xs,ωs) do x,ω
        abs(pf(x,ω) - sf(x,ω))
    end

    @test maximum(errs) < 2.0 * eps(Float64)

    # test constructors for 3D acoustics checks
    a3_host = Acoustic(3, ρ = 1.3, c = 1.5 + 0.0im)

    direction = rand(3)
    position = rand(3)
    amplitude = rand() + 0.1

    # check expansion in regular spherical waves
    psource = plane_source(a3_host; direction = direction, position =  position, amplitude = amplitude)
    pointsource = point_source(a3_host, position, amplitude)


    # Check that the field converges to its regular basis expansion around centre
    k = real(ω / psource.medium.c)
    centre = position + (4.0pi / k) .* rand(3)
    x = centre + (0.1 / k) .* rand(3)

    s_expand = source_expand(psource, x; basis_order = 4)

    # test if the source is equal to its series expansion in a regular basis.
    @test norm(field(psource,x,ω) - s_expand(x,ω)) < 1e-10 * abs(field(psource,x,ω))

    # create source through a regular expansion of spherical waves

    order = 3
    len = basisorder_to_basislength(typeof(a3_host),order)

    position = rand(3)

    source = regular_spherical_source(a3_host, rand(len) ./ (1:len);
        amplitude = rand(Complex{Float64}),
        position = position
    );

    k = real(ω / source.medium.c)
    centre = position + (4.0pi / k) .* rand(3)
    x = centre + (0.1 / k) .* rand(3)

    s_expand = source_expand(source, x; basis_order = 4)

    @test norm(field(source,x,ω) - s_expand(x,ω)) < 1e-10 * abs(field(source,x,ω))

end
