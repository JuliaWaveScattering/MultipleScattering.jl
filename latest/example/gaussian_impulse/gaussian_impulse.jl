using MultipleScattering
using Plots
pyplot()

# A very sharp sinc function (dicrete delta)
w_arr = collect(linspace(0.0,1.0,1000))
x0 = 10.;
t_arr = 4.:0.01:22;
fs = reshape([exp(im*w*x0) for w in w_arr],length(w_arr),1)
fts = frequency_to_time(fs,w_arr, t_arr);

# A gaussian covering same frequency range which leads to a time pulse with height 1
a = 3./maximum(w_arr)^2;
gauss_impulse = get_gaussian_freq_impulse(maximum(w_arr))
ft2s = frequency_to_time(fs,w_arr, t_arr; impulse = gauss_impulse);

# true impulse in time
true_fts = [exp(-(t-x0)^2/(4a)) for t in t_arr]
plot(t_arr, [fts,true_fts,ft2s])
