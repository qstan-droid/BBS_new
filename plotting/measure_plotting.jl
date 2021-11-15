using DelimitedFiles
using Plots
import PyPlot
const plt = PyPlot

# import the data
parent_folder = "one mode stuff"
folder_names = [string(parent_folder, "\\data_ml_ave_0.01_same_err_N2\\")] # put folder names in here

gate_het = 1 .- readdlm(string(@__DIR__, "\\", folder_names[1], "gate_het.csv"), ',', Float64)
error_bars_het = readdlm(string(@__DIR__, "\\", folder_names[1], "gate_SE_het.csv"), ',', Float64)

gate_ahd = 1 .- readdlm(string(@__DIR__, "\\", folder_names[1], "gate_ahd.csv"), ',', Float64)
error_bars_ahd = readdlm(string(@__DIR__, "\\", folder_names[1], "gate_SE_ahd.csv"), ',', Float64)

gate_opt = 1 .- readdlm(string(@__DIR__, "\\", folder_names[1], "gate_opt.csv"), ',', Float64)
error_bars_opt = readdlm(string(@__DIR__, "\\", folder_names[1], "gate_SE_opt.csv"), ',', Float64)

# ranges to plot against
no_of_samples = 1000
N_ord = 2

K = 1:12
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
    
p = plot(n_ave, 
        [gate_het gate_ahd gate_opt],
        yerr = [error_bars_het error_bars_ahd error_bars_opt],
        color = [:red :blue :green],
        xlabel = "\$\\hat{n}\$",
        ylabel = "Average Infidelity",
        label = ["HET" "AHD" "OPT"],
        title = "nu = 0.01 on both qubits | N = 2 | ml ave decoder",
        yaxis = :log,
        ylim = (10^-3, 1),
        titlefontsize = 8)
#Plots.abline!(0, break_even, line = :dash, label = false)

#p = plot(n_ave, 
#        [gate_fid[1] gate_fid[2]],
#        yerr = [err_bars[1] err_bars[2]],
#        color = [:red :blue],
#        xlabel = "\$\\hat{n}\$",
#        ylabel = "Average Infidelity",
#        label = ["naive" "ml ave"],
#        title = "\$ \nu = 0.01\$ on both rails | 3-qubit rep code | N = 2",
#        yaxis = :log,
#        titlefontsize = 8)

savefig(p, string(@__DIR__, "\\important plots\\measure_compare_0.01_N2.png"))