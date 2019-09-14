# function plot{T}(mnts::StatisticalMoments{T}; ribbon=false, X = mnts.x_arr, use_moments = [1,2,4]) #, kws...)
#     # if ribbon == true
#       # plot_ribbon(mnts)
#       plot_ribbon(mnts; X=X, use_moments=use_moments)
#     # else
#     #   plot_graph(mnts; kws...)
#     # end
#     # println("well hello mister put me nose in other peoples code")
# end

# # woud have liked a user series plot:
# """
#     plot(x, y, seriestype = :moments; kws...)
#
# A failed attempt a creating a user defined series plot for the moments `y` with x-axis `x`.
# """
@recipe function plot_ribbon(results::AbstractVector{SimRes};
        field_apply = real,
        num_moments = 2,
        use_moments = 1:min(num_moments,length(results)),
        Y = statistical_moments(results, num_moments; field_apply=field_apply)[use_moments],
        x = results[1].ω,
        labs = ["mean" "std" "skew" "kurt"]) where {T,SimRes<:SimulationResult{T}}

    # must include mean
    use_moments = sort(union([1; use_moments]))
    if labs == [] labs = repeat([""], inner = length(y)) end
    colors =[:black, :green, :orange, :red]

    m = Y[1]
    Ys = map(2:length(Y)) do i
        Y_down = if iseven(use_moments[i])
            Y[i]
            # Y[i]/T(2)
        else
            map(y-> (y< zero(T)) ? -y : zero(T), Y[i])
        end
        Y_up = if iseven(use_moments[i])
            Y[i]
            # Y[i]/T(2)
        else
            map(y-> (y> zero(T)) ? y : zero(T), Y[i])
        end
        (Y_down, Y_up)
    end
    for i = 1:(length(Ys)-1)
        Ys[i+1] = (Ys[i+1][1] + Ys[i][1], Ys[i+1][2] + Ys[i][2])
    end

    for i in reverse(1:length(Ys))
      @series begin
        seriestype := :line
        ribbon := Ys[i]
        lab --> labs[use_moments[i+1]]
        c --> colors[use_moments[i+1]]
        (x, m[:])
      end
    end

    for i in 1:min(15,length(results))
      @series begin
        linealpha --> 0.2
        lab --> ""
        c --> :grey
        results[i]
      end
    end

    @series begin
      seriestype := :line
      c --> :black
      linewidth --> 2
      lab --> labs[1]
      (x, m[:])
    end
end
