@testset "Symmetries" begin

    rectangle = Box([[0.0,0.0], [2.0, 3.0]])
    circle = Sphere(rand(2), 2.0)
    h2Dazi = Halfspace([1.0,0.0], [0.0,0.0])

    shapes2D = [rectangle,circle,h2Dazi]

    a2 = Acoustic(0.1,0.1 + 0.0im,2)

    source_position = [0.0,1.0]
    s1 = point_source(a2, source_position)
    s2 = plane_source(a2; position = source_position, direction = [1.0,0.0])

    # NOTE: should use the function regular_spherical_source, but for 2D one of the translation matrices has not been implemented yet
    s3 = RegularSource{Acoustic{Float64,2},RadialSymmetry{2}}(a2, (x,ω) -> norm(x)+1.0im, (order,centre,ω) -> [1.0+1.0im for m = -order:order])

    sources2D = [s1,s2,s3]

    syms = [Symmetry(sh,s) for sh in shapes2D, s in sources2D]
    @test syms == AbstractSymmetry{2}[
        WithoutSymmetry{2}() WithoutSymmetry{2}() WithoutSymmetry{2}();
        WithoutSymmetry{2}() AzimuthalSymmetry{2}() RadialSymmetry{2}();
        WithoutSymmetry{2}() PlanarAzimuthalSymmetry{2}() AzimuthalSymmetry{2}()
    ]

    check_commute = [Symmetry(sh,s) == Symmetry(s,sh) for s in sources2D, sh in shapes2D]
    @test count(check_commute) == length(sources2D) * length(shapes2D)


    sphere = Sphere(rand(3), 2.0)
    p = Plate([1.0,1.0,1.0], 1.0)
    hazi = Halfspace([0.0,0.0,1.0])

    shapes3D = [p,sphere,hazi]

    a3 = Acoustic(3, ρ = 1.3, c = 1.5)

    psource = plane_source(a3; direction = [0.0,0.0,1.0])
    pointsource = point_source(a3, [0.0,0.0,0.0])
    radsource = regular_spherical_source(a3, [1.0+0.0im];
       symmetry = RadialSymmetry{3}()
    );

    sources3D = [pointsource,psource,radsource];

    syms = [Symmetry(sh,s) for sh in shapes3D, s in sources3D];
    @test syms == AbstractSymmetry{3}[
        WithoutSymmetry{3}() PlanarSymmetry{3}()  WithoutSymmetry{3}();
        WithoutSymmetry{3}() AzimuthalSymmetry{3}() RadialSymmetry{3}();
        WithoutSymmetry{3}() PlanarAzimuthalSymmetry{3}() AzimuthalSymmetry{3}()
    ]

    check_commute = [Symmetry(sh,s) == Symmetry(s,sh) for s in sources3D, sh in shapes3D]
    @test count(check_commute) == length(sources3D) * length(shapes3D)

end
