# function plot{T}(mnts::Moments{T}; ribbon=false, X = mnts.x_arr, use_moments = [1,2,4]) #, kws...)
#     # if ribbon == true
#       # plot_ribbon(mnts)
#       plot_ribbon(mnts; X=X, use_moments=use_moments)
#     # else
#     #   plot_graph(mnts; kws...)
#     # end
#     # println("well hello mister put me nose in other peoples code")
# end

@recipe function plot_ribbon{T}(mnts::Moments{T};
      use_moments = [1,2,4], X = mnts.x_arr)

    labs = ["mean" "std" "skew" "kurt"]
    use_moments = sort(union([1; use_moments]))
    if labs == [] labs = repeat([""], inner = length(mnts.moments)) end
    colors =[:black, :green, :orange, :red]

    Y = mnts.moments[use_moments]
    m = Y[1]
    Y_up = zero(T)*m
    Y_down = zero(T)*m
    Ys = Array{Tuple}(length(Y)-1)
    for i = 2:length(Y)
      begin
        if iseven(use_moments[i])
          Y_up = Y_up + Y[i]/(2*one(T))
          Y_down = Y_down + Y[i]/(2*one(T))
        else
          Y_up = Y_up + map(y-> (y> zero(T))? y : zero(T), Y[i])
          Y_down = Y_down + map(y-> (y< zero(T))? -y : zero(T), Y[i])
        end
        Ys[i-1] = (Y_down, Y_up)
      end
    end
    for i in reverse(1:length(Ys))
      @series begin
        ribbon := Ys[i]
        lab --> labs[use_moments[i+1]]
        c --> colors[use_moments[i+1]]
        (X, m)
      end
    end
    @series begin
      c --> :black
      linewidth --> 2
      lab --> labs[1]
      (X, m)
    end
end
