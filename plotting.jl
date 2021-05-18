using Plots
using DelimitedFiles

folder = "data_single_mode_error_(-1)"
het = "heterodyne"
opt = "opt_phase"

fid_het = readdlm(string("C:\\Users\\Edgar\\Desktop\\Compare With Old Code Plotting\\", folder, "\\", het,".csv"), ',', Float64)
fid_opt = readdlm(string("C:\\Users\\Edgar\\Desktop\\Compare With Old Code Plotting\\", folder, "\\", opt,".csv"), ',', Float64)

###### Plotting ######

# First we get the average 
code = "binomial"
N_ord = 3

if code == "binomial"
    K = 1:1:15
    n_ave = N_ord .* K ./ 2
elseif code == "cat"
    alpha = 0.1:0.1:3
end

# We also calculate errorbars
no_samples = 1000
fid_het_error = sqrt.(fid_het .* (1 .- fid_het) ./ no_samples)
fid_opt_error = sqrt.(fid_opt .* (1 .- fid_opt) ./ no_samples)

# we then calculate break-even
break_even = function(nu)
    1 - (1 + exp(-nu/2))^2/4
end

nu = 10^(-1)
row, col = size(fid_het)
break_even_plot = zeros(Float64, (row, col))
for i = 1:row
    break_even_plot[i] = break_even(nu)
end

break_even_error = zeros(Integer, (row, col))

##############################

println([(1 .- fid_het) (1 .- fid_opt)])
println(1 .- fid_het)

p = plot(n_ave,
            [(1 .- fid_het) (1 .- fid_opt)],
            yerror = [fid_het_error fid_opt_error],
            shape = [:circle :circle],
            linestyle = [:solid :solid],
            markersize = 2,
            color = [:blue :red],
            xlabel = "\$\\hat{n}\$",
            ylabel = "Average Infidelity",
            label = ["heterodyne" "opt_phase"],
            yaxis = :log,
            title = "loss rate = 10^{-1} | single mode | order = [3, 1] | K_{2} = 25",
            titlefontsize = 11,
            ylim = (10^(-2), 1))

plot!(n_ave, break_even_plot, linestyle = :dash, label = "break even")

display(p)
savefig(p, string("C:\\Users\\Edgar\\Desktop\\Compare With Old Code Plotting\\", folder, "\\plot.png"))