using DelimitedFiles
using Plots
import PyPlot
const plt = PyPlot

# import the data
parent_folder = "rep mode stuff"
folder_names = [string(parent_folder, "\\data_1_3_1_ml_ave_0.01_same_err_no_diff_N1_opt_more_samp"), string(parent_folder, "\\data_1_3_1_ml_ave_0.01_same_err_no_diff_N2_opt_more_samp")] # put folder names in here

gate_fid = []
err_bars = []

gate_fid_opt = []
err_bars_opt = []

# heterodyne or opt_phase
#het = "\\gate_het.csv"
het = "\\gate_opt.csv"
#het = "\\gate_ahd.csv"

#het_SE = "\\gate_SE_het.csv"
het_SE = "\\gate_SE_opt.csv"
#het_SE = "\\gate_SE_ahd.csv"

for i = 1:length(folder_names)
    gate_sum_import = 1 .- readdlm(string(@__DIR__, "\\", folder_names[i], het), ',', Float64)
    error_bars_import = readdlm(string(@__DIR__, "\\", folder_names[i], het_SE), ',', Float64)

    #gate_sum_import_opt = 1 .- readdlm(string(@__DIR__, "\\", folder_names[i], opt), ',', Float64)
    #error_bars_import_opt = readdlm(string(@__DIR__, "\\", folder_names[i], opt_SE), ',', Float64)

    push!(gate_fid, gate_sum_import)
    push!(err_bars, error_bars_import)

    #push!(gate_fid_opt, gate_sum_import_opt)
    #push!(err_bars_opt, error_bars_import_opt)
end

# ranges to plot against
no_of_samples = 5000
N_ord = 1

K = 1:13
n_ave = N_ord .* K ./ 2

##### Calculating Break Even #####
nu = 0.01
break_even = 1 - ((1 + exp(-nu/2))^2 /2 + 1)/3

#################################
#p = plot(n_ave, 
#        [gate_fid_opt[1] gate_fid_opt[2]],
#        yerr = [err_bars_opt[1] err_bars_opt[2]],
#        color = [:red :blue],
#        xlabel = "\$\\hat{n}\$",
#        ylabel = "Average Infidelity",
#        label = ["no rep" "rep = 3"],
#        title = "\$ \nu = 0.01\$ on both rails | optimal phase measurements | average ML decoder | N = 1",
#        yaxis = :log,
#        titlefontsize = 8)
    
p = plot(K, 
        [gate_fid[1][K] gate_fid[2][K]],
        yerr = [err_bars[1][K] err_bars[2][K]],
        color = [:red :blue],
        xlabel = "K",
        ylabel = "Average Infidelity",
        label = ["N = 1" "N = 2"],
        title = string("nu = ", nu, " both blocks | opt measurements | ml ave decoder | m = 3, n = 1, rep = 1 | N = 5000"),
        yaxis = :log,
        #ylim = (10^-1, 1),
        titlefontsize = 8)
Plots.abline!(0, break_even, line = :dash, label = false)

#p = plot(n_ave, 
#        [gate_fid[1] gate_fid[2]],
#        yerr = [err_bars[1] err_bars[2]],
#        color = [:red :blue],
#        xlabel = "\$\\hat{n}\$",
#        ylabel = "Average Infidelity",
#        label = ["naive" "ml ave"],
#        title = "\$ \nu = 0.01\$ on both rails | 3-qubit rep code | N = 1",
#        yaxis = :log,
#        titlefontsize = 8)\
plot(p)
savefig(p, string(@__DIR__, "\\important plots\\rep_order_comparem=3_0.01_opt_more_samp.png"))