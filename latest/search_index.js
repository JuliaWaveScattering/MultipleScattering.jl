var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#MultipleScattering:-wave-scattering-from-multiple-particles-in-Julia-1",
    "page": "Home",
    "title": "MultipleScattering: wave scattering from multiple particles in Julia",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Install Julia v0.6.0 or later, if you haven\'t already.Pkg.clone(\"https://github.com/jondea/MultipleScattering.jl.git\")\nusing MultipleScattering"
},

{
    "location": "index.html#Contents-1",
    "page": "Home",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"base.md\", \"acoustics.md\",\"random.md\"]\nDepth = 3"
},

{
    "location": "index.html#Index-1",
    "page": "Home",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "base.html#",
    "page": "Base",
    "title": "Base",
    "category": "page",
    "text": ""
},

{
    "location": "base.html#Base-1",
    "page": "Base",
    "title": "Base",
    "category": "section",
    "text": ""
},

{
    "location": "base.html#MultipleScattering.Shape",
    "page": "Base",
    "title": "MultipleScattering.Shape",
    "category": "type",
    "text": "Abstract idea which defines the external boundary of object.\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.origin",
    "page": "Base",
    "title": "MultipleScattering.origin",
    "category": "function",
    "text": "origin(shape::Shape)::SVector\n\nOrigin of shape, typically the center\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.iscongruent",
    "page": "Base",
    "title": "MultipleScattering.iscongruent",
    "category": "function",
    "text": "iscongruent(p1::Shape, p2::Shape)::Bool\n\nTrue if shapes are the same but in different positions (origins), standard mathematical definition.\n\n\n\niscongruent(p1::AbstractParticle, p2::AbstractParticle)::Bool\n\nReturns true if medium and shape of particles are the same, ignoring origin, false otherwise.\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.congruent",
    "page": "Base",
    "title": "MultipleScattering.congruent",
    "category": "function",
    "text": "congruent(s::Shape, x)::Shape\n\nCreate shape congruent to s but with origin at x\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.bounding_rectangle",
    "page": "Base",
    "title": "MultipleScattering.bounding_rectangle",
    "category": "function",
    "text": "Returns rectangle which completely encloses the shapes\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.boundary_functions",
    "page": "Base",
    "title": "MultipleScattering.boundary_functions",
    "category": "function",
    "text": "volume(shape::Shape)::NTuple{Function,Dim)\n\nReturns Tuple of Dim Functions which define outer boundary of shape when given boundary coordinate t∈[0,1]\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.name",
    "page": "Base",
    "title": "MultipleScattering.name",
    "category": "function",
    "text": "name(shape::Shape)::String\n\nName of a shape\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.outer_radius",
    "page": "Base",
    "title": "MultipleScattering.outer_radius",
    "category": "function",
    "text": "outer_radius(shape::Shape{T})::T\n\nThe radius of a circle which completely contains the shape\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.volume",
    "page": "Base",
    "title": "MultipleScattering.volume",
    "category": "function",
    "text": "volume(shape::Shape{T})::T\n\nVolume of a shape\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.Circle",
    "page": "Base",
    "title": "MultipleScattering.Circle",
    "category": "type",
    "text": "Circle(origin, radius)\n\n2D Shape where boundary is a fixed distance from the origin\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.Rectangle",
    "page": "Base",
    "title": "MultipleScattering.Rectangle",
    "category": "type",
    "text": "Rectangle(origin::AbstractVector{T}, width::T, Height::T)\nRectangle(bottomleft::AbstractVector{T}, topright::AbstractVector{T})\n\n2D Shape with axis aligned sides, defined by width, height and origin (at the center).\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.TimeOfFlightFromPoint",
    "page": "Base",
    "title": "MultipleScattering.TimeOfFlightFromPoint",
    "category": "type",
    "text": "A shape where anything inside could cause a disturbance at the listener position from a point source wavefront starting at the listener. Also everything inside has a positive x coordinate. It is equivalent to a segment of a circle.\n\nMore precisely, if the listener is at (l_x,l_y) then the interior of the shape is defined as sqrt((x-l_x)^2+(y-l_y)^2)<time and x>0\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.TimeOfFlight",
    "page": "Base",
    "title": "MultipleScattering.TimeOfFlight",
    "category": "type",
    "text": "A shape where anything inside could cause a disturbance at the listener position from a planar wavefront parallel to the y axis starting at the listener. Also everything inside has a positive x coordinate.\n\nMore precisely, if the listener is at (l_x,l_y) then the interior of the shape is defined as x-l_x+sqrt((x-l_x)^2+(y-l_y)^2)<time and x>0\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.Sphere",
    "page": "Base",
    "title": "MultipleScattering.Sphere",
    "category": "type",
    "text": "Shape where boundary is a fixed distance from the origin\n\n\n\n"
},

{
    "location": "base.html#Shapes-1",
    "page": "Base",
    "title": "Shapes",
    "category": "section",
    "text": "Shape types and functions.MultipleScattering.Shape\nMultipleScattering.origin\nMultipleScattering.iscongruent\nMultipleScattering.congruent\nMultipleScattering.bounding_rectangle\nMultipleScattering.boundary_functions\nMultipleScattering.name\nMultipleScattering.outer_radius\nMultipleScattering.volume\nMultipleScattering.Circle\nMultipleScattering.Rectangle\nMultipleScattering.TimeOfFlightFromPoint\nMultipleScattering.TimeOfFlight\nMultipleScattering.Sphere"
},

{
    "location": "base.html#MultipleScattering.PhysicalProperties",
    "page": "Base",
    "title": "MultipleScattering.PhysicalProperties",
    "category": "type",
    "text": "Holds information about the physical properties of the medium, the dimension of the field and the number of dimensions it is a function of.\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.field_dim",
    "page": "Base",
    "title": "MultipleScattering.field_dim",
    "category": "function",
    "text": "Extract the dimension of the field of this physical property\n\n\n\nExtract the dimension of the field of this type of physical property\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.dim",
    "page": "Base",
    "title": "MultipleScattering.dim",
    "category": "function",
    "text": "Extract the dimension of the space that this physical property lives in\n\n\n\nExtract the dimension of the space that this type of physical property lives in\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.basis_function",
    "page": "Base",
    "title": "MultipleScattering.basis_function",
    "category": "function",
    "text": "Basis functions in a specific dimension for a specific physics type.\n\n\n\nBasis function when inside a particle. Assumes particle is a circle, which approximately works for all shapes.\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.internal_field",
    "page": "Base",
    "title": "MultipleScattering.internal_field",
    "category": "function",
    "text": "the field inside an AbstractParticle a some given point x.\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.boundary_data",
    "page": "Base",
    "title": "MultipleScattering.boundary_data",
    "category": "function",
    "text": "A tuples of vectors of the field close to the boundary of the shape. The field is calculated from sim::FrequencySimulation, but the PhysicalProperties inside and outside of the shape are assumed to be given by inside_medium and outside_medium.\n\n\n\n"
},

{
    "location": "base.html#Physical-properties-1",
    "page": "Base",
    "title": "Physical properties",
    "category": "section",
    "text": "Physical properties types and functions.MultipleScattering.PhysicalProperties\nMultipleScattering.field_dim\nMultipleScattering.dim\nMultipleScattering.basis_function\nMultipleScattering.internal_field\nMultipleScattering.boundary_data"
},

{
    "location": "base.html#Particles-1",
    "page": "Base",
    "title": "Particles",
    "category": "section",
    "text": "Particle types and functions.MultipleScattering.AbstractParticle\nMultipleScattering.Particle\nMultipleScattering.CapsuleParticle\nMultipleScattering.iscongruent"
},

{
    "location": "base.html#MultipleScattering.Source",
    "page": "Base",
    "title": "MultipleScattering.Source",
    "category": "type",
    "text": "Source{P,T}(field::Function,coef::Function) where P<:PhysicalProperties{T,Dim,FieldDim}\n\nA struct which describes the source that drives/forces the whole system. It is also described as an incident wave.\n\nx = [1.0,0.0]\nω = 1.0\nSource.field(x,ω)\n\nshould return the wave field at position \'x\' and angular frequency \'ω\'\n\n\n\n"
},

{
    "location": "base.html#MultipleScattering.self_test",
    "page": "Base",
    "title": "MultipleScattering.self_test",
    "category": "function",
    "text": "Check that the source functions return the correct types\n\n\n\nCheck that the discrete impulse vectors are the right sizes\n\n\n\nCheck that the continuous impulse functions return the correct types\n\n\n\n"
},

{
    "location": "base.html#Source-1",
    "page": "Base",
    "title": "Source",
    "category": "section",
    "text": "Source types and functions.MultipleScattering.Source\nMultipleScattering.self_test"
},

{
    "location": "base.html#Simulation-1",
    "page": "Base",
    "title": "Simulation",
    "category": "section",
    "text": "Simulation types and functions.MultipleScattering.Simulation\nMultipleScattering.FrequencySimulation\nMultipleScattering.run\nMultipleScattering.forcing\nMultipleScattering.basis_coefficients\nMultipleScattering.field\nMultipleScattering.scattering_matrix\nMultipleScattering.t_matrix\nMultipleScattering.get_t_matrices"
},

{
    "location": "acoustics.html#",
    "page": "Acoustic",
    "title": "Acoustic",
    "category": "page",
    "text": ""
},

{
    "location": "acoustics.html#Acoustic-1",
    "page": "Acoustic",
    "title": "Acoustic",
    "category": "section",
    "text": "Acoustic type and functions.MultipleScattering.Acoustic\nMultipleScattering.impedance\nMultipleScattering.sound_hard\nMultipleScattering.hard\nMultipleScattering.rigid\nMultipleScattering.zero_neumann\nMultipleScattering.sound_soft\nMultipleScattering.soft\nMultipleScattering.pressure_release\nMultipleScattering.zero_dirichlet\nMultipleScattering.plane_source\nMultipleScattering.point_source"
},

{
    "location": "random.html#",
    "page": "Random",
    "title": "Random",
    "category": "page",
    "text": ""
},

{
    "location": "random.html#Random-1",
    "page": "Random",
    "title": "Random",
    "category": "section",
    "text": ""
},

{
    "location": "random.html#MultipleScattering.random_particles",
    "page": "Random",
    "title": "MultipleScattering.random_particles",
    "category": "function",
    "text": "random_particles(particle_medium, particle_shape, box_shape, volume_fraction::Number;\n    seed=Base.Random.make_seed())\nrandom_particles(particle_medium, particle_shape, box_shape, N::Integer;\n    seed=Base.Random.make_seed())\n\nGenerate N random particles that fit inside box_shape (or fill with volume_fraction)\n\nSpecify seed to make output deterministic. Algorithm places particles unifomly randomly inside the bounding rectangle of box_shape and discards particle if it overlaps (based on outer radius) or does not lies completely in box box.\n\n\n\n"
},

{
    "location": "random.html#Random-particles-1",
    "page": "Random",
    "title": "Random particles",
    "category": "section",
    "text": "MultipleScattering.random_particles"
},

{
    "location": "random.html#MultipleScattering.calculate_moments",
    "page": "Random",
    "title": "MultipleScattering.calculate_moments",
    "category": "function",
    "text": "calculate_moments(results, n; applytofield=real)::Vector{Matrix}\n\nCalculate moments up to n of results at each position and wavenumber/time, after applying applytofield.\n\n\n\n"
},

{
    "location": "random.html#Moments-1",
    "page": "Random",
    "title": "Moments",
    "category": "section",
    "text": "MultipleScattering.calculate_moments"
},

{
    "location": "example/box_size/README.html#",
    "page": "Box size",
    "title": "Box size",
    "category": "page",
    "text": ""
},

{
    "location": "example/box_size/README.html#Box-size-1",
    "page": "Box size",
    "title": "Box size",
    "category": "section",
    "text": "If we are only interested in the response at a specific location for a certain time interval, we need only simulate particles which are this distance (or time) away from the listener."
},

{
    "location": "example/hankel_convergence/README.html#",
    "page": "Convergence when increasing the number of Hankel functions",
    "title": "Convergence when increasing the number of Hankel functions",
    "category": "page",
    "text": ""
},

{
    "location": "example/hankel_convergence/README.html#Convergence-when-increasing-the-number-of-Hankel-functions-1",
    "page": "Convergence when increasing the number of Hankel functions",
    "title": "Convergence when increasing the number of Hankel functions",
    "category": "section",
    "text": "The code convergence.jl tests how fast does the scattered wave (in frequency) converge when increasing the number of Hankel functions. To describe the scattered wave from each particle we use a series of Hankel functions (of the first kind).include(\"convergence.jl\")\nsimulations = hankel_order_convergence()\nplot_hankel_order_convergence(simulations)\n(Image: Plot lens shape and response in time)In the figures above m is the maximum order of the Hankel functions. The top left figure shows the configuration of particles considered.  "
},

{
    "location": "example/intro/README.html#",
    "page": "Go to introductory example.",
    "title": "Go to introductory example.",
    "category": "page",
    "text": ""
},

{
    "location": "example/intro/README.html#[Go-to-introductory-example](https://github.com/jondea/MultipleScattering.jl).-1",
    "page": "Go to introductory example.",
    "title": "Go to introductory example.",
    "category": "section",
    "text": ""
},

{
    "location": "example/lens/README.html#",
    "page": "Lens",
    "title": "Lens",
    "category": "page",
    "text": ""
},

{
    "location": "example/lens/README.html#Lens-1",
    "page": "Lens",
    "title": "Lens",
    "category": "section",
    "text": "the code lens.jl arranges particles into the shape of lens. The lens is shaped so that the incident plane wave is completely focused into one point: the listener position.include(\"lens.jl\")\nplot_lens()(Image: Plot lens shape and response in time)The distance of the lens is chosen so that the peak signal should arrive around time 34. Note, as we use a random seed to position the particles, the above figure many vary slightly."
},

{
    "location": "example/moments/README.html#",
    "page": "StatisticalMoments",
    "title": "StatisticalMoments",
    "category": "page",
    "text": ""
},

{
    "location": "example/moments/README.html#StatisticalMoments-1",
    "page": "StatisticalMoments",
    "title": "StatisticalMoments",
    "category": "section",
    "text": "Here we are going to simulate the scattered wave for many different configurations of particles. We can then take the average and standard deviation (the moments) of the scattered wave. In statistical mechanics this process is called ensemble average."
},

{
    "location": "example/moments/README.html#Choose-the-type-of-particles-1",
    "page": "StatisticalMoments",
    "title": "Choose the type of particles",
    "category": "section",
    "text": "using MultipleScattering\n\nvolfrac = 0.01\nradius = 1.0\nnum_particles = 10\n\n# region to place the particles\nshape = Rectangle(volfrac, radius, num_particles)To see the region where the particle will be placed, and the receiver position:using Plots\npyplot()\nlistener = [-10.0, 0.0]\nplot(shape);\nscatter!([listener[1]],[listener[2]]);\nplot_shape = annotate!([(listener_position[1], listener_position[2] -2., \"Receiver\")])(Image: Plot of shape and receiver)"
},

{
    "location": "example/moments/README.html#Calculate-the-moments-of-the-scattered-wave-1",
    "page": "StatisticalMoments",
    "title": "Calculate the moments of the scattered wave",
    "category": "section",
    "text": "The code below chooses a random (uniform distribution) configuration of particles inside shape and calculates the received signal at listener for wavenumbers k_arr,k_arr = collect(linspace(0.01,1.0,100))\nsimulation = FrequencySimulation(volfrac,radius,k_arr; shape=shape, listener_positions = listener, seed = 1)\nplot(simulation)(Image: Plot of response against wavenumber)To see the position of the chosen particles:plot(plot_shape)\nplot!.(simulation.particles);\nplot!()(Image: Plot of the position of the particles)Now we will do simulations for particles placed in many different configurations and take the moments:simulations = [\n    FrequencySimulation(volfrac,radius,k_arr; shape=shape, listener_positions = listener, seed = i)\nfor i = 1:20]\nreal_moments = StatisticalMoments(simulations; response_apply=real) # moments of the real part\nplot(real_moments);\nplot!(xlabel=\"wavenumbers\", title=\"Moments of the real part\")(Image: Moments of the real part the scattered waves)"
},

{
    "location": "example/moments/README.html#Calculate-the-moments-of-the-scattered-wave-in-time-1",
    "page": "StatisticalMoments",
    "title": "Calculate the moments of the scattered wave in time",
    "category": "section",
    "text": "time_simulations = TimeSimulation.(simulations)\ntime_simulations[1].time_arr # the time_arr chosen will be based on the discrete Fourier transform of simulations[1].k_arr\nreal_time_moments = StatisticalMoments(time_simulations; response_apply=real) # moments of the real part\nplot(real_time_moments,xlims=(0,300));\nplot!(xlabel=\"time\", title=\"Moments of the real part of the time wave\")(Image: Moments of the real part the scattered waves in time)"
},

{
    "location": "example/moments/README.html#References-1",
    "page": "StatisticalMoments",
    "title": "References",
    "category": "section",
    "text": "A. L. Gower, R. M. Gower, J. Deakin, W. J. Parnell, I. D. Abrahams, Learning about random media from near-surface backscattering: using machine learning to measure particle size and concentration, arXiv preprint, (2018)1801.05490"
},

{
    "location": "example/near_surface_backscattering/README.html#",
    "page": "Near-surface backscattering",
    "title": "Near-surface backscattering",
    "category": "page",
    "text": ""
},

{
    "location": "example/near_surface_backscattering/README.html#Near-surface-backscattering-1",
    "page": "Near-surface backscattering",
    "title": "Near-surface backscattering",
    "category": "section",
    "text": "Near-surface backscattering is a method of accurately calculating the backscattering from an infinite halfspace. First, let us see why it is difficult to approximate the scattering from a halfspace filled with particles. That is, let us find out how many particles are required before the backscattering converges."
},

{
    "location": "example/near_surface_backscattering/README.html#Generate-a-large-material-filled-with-particles.-1",
    "page": "Near-surface backscattering",
    "title": "Generate a large material filled with particles.",
    "category": "section",
    "text": "using MultipleScattering\nusing Plots\npyplot(linewidth=2)\n\nradius = 0.8\nvolfrac = 0.10\nmax_width = 60.\n\nbottomleft = [0.,-max_width]\ntopright = [max_width,max_width]\n\nshape = Rectangle(bottomleft,topright)\nparticles = random_particles(volfrac, radius, shape; c=1.0+0.0im, ρ=0.0)We will measure the backscattering at listener_position:listener_position = [-10.,0.]\nscatter([listener_position[1]],[listener_position[2]]);\nannotate!([(listener_position[1], listener_position[2] -max_width/10., \"Receiver\")])\nplot!.(particles);\nplot!(shape)(Image: The largest quantity of particles used)"
},

{
    "location": "example/near_surface_backscattering/README.html#Calculate-backscattering-for-different-quantity-of-particles-1",
    "page": "Near-surface backscattering",
    "title": "Calculate backscattering for different quantity of particles",
    "category": "section",
    "text": "We will shave off particles on the right of this group of particles (above), and then calculate the resulting backscattered waves.widths = 10.:10.:max_width # choose the width of the region filled with particles\nk_arr = collect(0.01:0.01:1.) # choose the wavenumbers of the incident wave\n\nsimulations = map(widths) do w # this is a for loop over the array widths\n    shape.topright[1] = w # choose a material with a smaller width\n    ps = filter(p -> p⊆shape, particles) # shave off particles\n    FrequencySimulation(ps, k_arr) # calculate backscattering\nend\n\nbackscattered_waves = [s.response for s in simulations]\nnum_particles = [length(s.particles) for s in simulations]\n\nM = length(backscattered_waves)\nbM = backscattered_waves[M] # backscattering from largest material\ndifferences = [norm(b - bM) for b in backscattered_waves[1:(M-1)]]./norm(bM)\n\nplot_converge = plot(num_particles[1:(M-1)], differences, xlabel = \"number of particles\", ylabel =\"error %\", label=\"frequency convergence\")(Image: The convergence of the response in frequency, when increasing the number of particles)The graph shows the rate of convergence, that is, how much the backscattering changes when including more particles (making the material deeper). The graph has not clearly converged, so we can only conclude that more than 400 particles are needed to accurately approximate the backscattering from an infinite halfspace. We can accelerate this convergence by considering backscattering in time."
},

{
    "location": "example/near_surface_backscattering/README.html#Calculate-backscattering-in-time-1",
    "page": "Near-surface backscattering",
    "title": "Calculate backscattering in time",
    "category": "section",
    "text": "time_simulations = TimeSimulation.(simulations)\n\ntimes = 2*(widths[1:5] .- listener_position[1]) # time if takes for an incident plane wave to reach the furthest particles and then return to the receiver\n\nplot()\nfor i in [1,3,6,9,12,13]\n    plot!(time_simulations[i],label=\"$(num_particles[i]) particles\"\n        , xlims=(0,maximum(times)+10.), ylims=(-0.6,0.3)\n        , xticks = [0.; 30.; times]\n    )\nend\ngui()(Image: The responses in time for different quantity of particles)We see that the responses in time diverge from each other more and more as time goes by. Meaning that if we only calculate the response for a short amount of time 34, then the convergence will be accelerated.time_arr = 0.:pi:34.2\ntime_simulations = [TimeSimulation(s;time_arr=time_arr) for s in simulations]\n\nbackscattered_waves = [s.response for s in time_simulations]\nbM = backscattered_waves[M] # backscattering from largest material\ndifferences = [norm(b - bM) for b in backscattered_waves[1:(M-1)]]./norm(bM)\nplot(plot_converge)\nplot!(num_particles[1:(M-1)], differences, xlabel = \"number of particles\", ylabel =\"error %\", label=\"time convergence\")(Image: Compare converges for responses in time and responses in frequency)The convergence of the time response, for time 0<t<34, is much faster. In fact, less than 100 particles are needed to accurately approximate the backscattering from an infinite halfspace. The reason we don\'t show these as log plots is because there is a small constant error (about 0.01%) due to the discrete Fourier transform. This error is caused by the Gibbs phenomena and by assuming the backscattering is periodic (when it is not). Both these errors are well understood and can be controlled."
},

{
    "location": "example/near_surface_backscattering/README.html#Calculate-backscattering-only-from-near-surface-particles-1",
    "page": "Near-surface backscattering",
    "title": "Calculate backscattering only from near-surface particles",
    "category": "section",
    "text": "This last step is about efficiency. We want to only include particle which contribute to the backscattering for short time intervals. To do this we created a region called TimeOfFlight(listener,time), where every particle in this shape takes less than time for their first scattered wave (due to an incident plane wave) to return to the listener.  More precisely, if listener = (lx,ly), then every point (x,y) inside this shape satisfies: x-lx+((x-lx)^2+(y-ly)^2)^(1/2)<time and x>0.For example, look at the largest quantity of particle we usedlistener_position = [-10.,0.]\nshape = TimeOfFlight(listener_position,80.0)\nscatter([listener_position[1]],[listener_position[2]]);\nannotate!([(listener_position[1], listener_position[2] -max_width/10., \"Receiver\")])\nplot!.(particles);\nplot!(shape)(Image: Shows the particles in the shape TimeOfFlight)For time 0<t<80 the backscattering from these particles is the same as an infinite halfspace filled with particles. To achieve this result we need only the particles inside the shape TimeOfFlight (region with the red outline). The particles outside this shape were unnecessary. To see this inaction:times = 40.:15.:80.\nnear_surface_simulations = map(times) do t\n    shape = TimeOfFlight(listener_position,t) # choose a material with particles only in the near surface region\n    ps = filter(p -> p⊆shape, particles) # shave off particles\n    FrequencySimulation(ps, k_arr; shape=shape) # calculate backscattering\nend\n\ntime_near_simulations = TimeSimulation.(near_surface_simulations)\n\nplot()\nfor i in 1:length(times)\n    plot!(time_near_simulations[i],label=\"time of flight $(times[i])\"\n        , xlims=(0,maximum(times)+10.), ylims=(-0.6,0.3)\n        , xticks = [0.; times], title=\"\"\n    )\nend\ngui()(Image: Response from particles in the shapes TimeOfFlight)Note the incident pulse has a thickness of about 10 in time, which is why the time of flight 40 diverges from the other curves slightly before time 40, and likewise for the other curves."
},

{
    "location": "example/random_particles/README.html#",
    "page": "Simple random particles example",
    "title": "Simple random particles example",
    "category": "page",
    "text": ""
},

{
    "location": "example/random_particles/README.html#Simple-random-particles-example-1",
    "page": "Simple random particles example",
    "title": "Simple random particles example",
    "category": "section",
    "text": "If it isn\'t installed, clone it from githubtry using MultipleScattering\ncatch\n    Pkg.clone(\"https://github.com/jondea/MultipleScattering.jl.git\")\nend\n\nusing MultipleScatteringDefine the fraction of the volume the particles will take up, their radius and which wavenumbers (k) to evaluate atvolfrac = 0.01\nradius = 1.0\nk_arr = collect(linspace(0.01,1.0,100))\nsimulation = FrequencySimulation(volfrac,radius,k_arr)We use the Plots package to plot both the response at the listener position and the whole field for a specific wavenumber (k=0.8)using Plots\nplot(\n    plot(simulation),\n    plot(simulation,0.8;res=100,drawshape=true)\n)(Image: Plot of response against wavenumber)(Image: Plot real part of acoustic field)"
},

{
    "location": "example/random_particles/README.html#Things-to-try-1",
    "page": "Simple random particles example",
    "title": "Things to try",
    "category": "section",
    "text": "Try changing the volume fraction, particle radius and k values we evaluate"
},

{
    "location": "example/random_particles_in_circle/README.html#",
    "page": "Random particles in a circle",
    "title": "Random particles in a circle",
    "category": "page",
    "text": ""
},

{
    "location": "example/random_particles_in_circle/README.html#Random-particles-in-a-circle-1",
    "page": "Random particles in a circle",
    "title": "Random particles in a circle",
    "category": "section",
    "text": "The code particles_in_circle.jl compares the scattered wave from one big circle, with the scattered wave from a circle filled with small particles.using MultipleScattering\n\nk_arr = collect(linspace(0.1,1.0,10))\n\n# You can also pick your own shape, an generate random particles inside it\n# with a certain radius ands volume fraction\nradius = 0.3\nvolfrac = 0.45\ncentre = [0.,0.]\nbig_radius = 3.0\n\ncircle = Circle(big_radius,centre)\ncircle_simulation = FrequencySimulation(volfrac,radius,k_arr;shape=circle)The particles chosen are impenetrable, i.e. the wave is 100\\% reflected. So this circle filled with scatterers should act like one big particle.big_particle = Particle(centre,big_radius)\nbig_particle_simulation = FrequencySimulation([big_particle], k_arr; hankel_order=15)\n\nusing Plots\nplot(\n    plot(circle_simulation,0.5; drawshape=true, drawlisteners =false),\n    plot(big_particle_simulation,0.5; drawlisteners =false),\n    layout = (2,1)\n)(Image: The field comparison)If we compare the response measured at the listener [-10., 0.], they should be very similar:plot(circle_simulation)\nplot!(big_particle_simulation,title=\"Compare scattered wave from one big particle, \\n and a circle filled with small particles\")(Image: The response comparison)"
},

{
    "location": "example/shapes/README.html#",
    "page": "Shapes",
    "title": "Shapes",
    "category": "page",
    "text": ""
},

{
    "location": "example/shapes/README.html#Shapes-1",
    "page": "Shapes",
    "title": "Shapes",
    "category": "section",
    "text": ""
},

{
    "location": "example/shapes/README.html#Existing-shapes-1",
    "page": "Shapes",
    "title": "Existing shapes",
    "category": "section",
    "text": "The package provides 3 built in basic shapes to put your random particles in, you can plot them using:using MultipleScattering\n\nrectangle = Rectangle([0.0,0.0],[2.0,3.0])\nplot(rectangle)\n\ncircle = Circle(3.0,[-1.0,2.0])\nplot!(circle)\n\ntimeofflight = TimeOfFlight([-1.0,0.0],3.0)\nplot!(timeofflight)Time of flight is a shape which contains shapes from a half space which take at most t time to reach from the listener."
},

{
    "location": "example/shapes/README.html#New-shape-1",
    "page": "Shapes",
    "title": "New shape",
    "category": "section",
    "text": "If you are feeling very adventurous, you can define your own shape First you must import the package in order to add to existing functionsimport MultipleScattering\n\ntype MyShape <: MultipleScattering.Shape\nendTo describe the characteristics and behaviour of the function you must define the following functions:MultipleScattering.volume(shape::MyShape) = 0.0\n\nMultipleScattering.name(shape::MyShape) = \"MyShape\"\n\nMultipleScattering.bounding_box(shape::MyShape) = MultipleScattering.Rectangle()When you have this, you can make use of your shape to generate particles in it"
},

{
    "location": "example/time_response_single_particle/README.html#",
    "page": "Time response from single particle",
    "title": "Time response from single particle",
    "category": "page",
    "text": ""
},

{
    "location": "example/time_response_single_particle/README.html#Time-response-from-single-particle-1",
    "page": "Time response from single particle",
    "title": "Time response from single particle",
    "category": "section",
    "text": ""
},

{
    "location": "example/two_particles/README.html#",
    "page": "Two particles",
    "title": "Two particles",
    "category": "page",
    "text": ""
},

{
    "location": "example/two_particles/README.html#Two-particles-1",
    "page": "Two particles",
    "title": "Two particles",
    "category": "section",
    "text": "Define two particles with the first centred at [1.,-2.], with radius 1.0, sound speed 2.0 and density 10.0using MultipleScattering\nusing Plots\npyplot()\n\np1 = Particle([1.,-4.], 1.0; c = 20.0+0.0im, ρ = 10.)\np2 = Particle([3.,3.],  3.0; c = 1.0+0.0im, ρ = 0.1)\nparticles = [p1,p2]Specify the angular frequency of the incident wave and calculate the responsew_arr = collect(0.1:0.01:1.)\nsimulation = FrequencySimulation(particles, w_arr)\nplot(simulation)(Image: Plot against frequency)The above used an incident plane with the default reciever/listener position and incident plane wave directionsimulation.listener_positions\nsimulation.source_directionto change these defaults usesimulation = FrequencySimulation(particles, w_arr;\n    listener_positions = [-10.,-10.],\n    source_direction=[1.,1.])then plot the response around the particles and receiverw = 3.2\nplot(simulation,w; res=80, resp_fnc=abs)(Image: Plot absolute value of wave field)the green circle in the plot is the receiver position. Looking at the region between the particles we see the complicated results of multiple scatttering."
},

]}
