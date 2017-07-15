# """
# Submodule to plot models built with MultipleScattering. Is a submodule due to
# size of package Plots, which is a requirement
# """
module Plot

export plot_particles, plot_shape, plot_listeners, plot_domain,
       build_field_model, plot_field_model, plot_field,
       plot_model,
       plot_moments

using MultipleScattering
import Plots: plot, plot!, pyplot, cgrad
import PlotUtils: ColorGradient
import ColorTypes: RGBA

include("domain.jl")
include("model.jl")
include("field.jl")
include("moments.jl")


end
