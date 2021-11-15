using DelimitedFiles

# which folder
parent_folder = "rep mode stuff"
folder = string(parent_folder, "\\data_1_3_1_ml_ave_0.01_same_err_no_diff_N2_opt_more_samp")

# Collecting all the files 
no_of_files = 10
test = readdlm(string(@__DIR__, "\\", folder, "\\2.csv"), ',', Float64)
l = length(test)

fidelity = zeros(l)
gate_fid = zeros(l)

# sum all fidelities
for i = 1:no_of_files
    for j = 1:l
        fidelity[j] = fidelity[j] + readdlm(string( @__DIR__, "\\", folder, "\\", i, ".csv"), ',', Float64)[j]
        gate_fid[j] = gate_fid[j] + readdlm(string( @__DIR__, "\\", folder, "\\gate_", i, ".csv"), ',', Float64)[j]
    end
end

# average all fidelities
for i = 1:l
    fidelity[i] = fidelity[i]/no_of_files 
    gate_fid[i] = gate_fid[i]/no_of_files 
end

writedlm(string(@__DIR__, "\\", folder, "\\opt.csv"), fidelity, ',')
writedlm(string(@__DIR__, "\\", folder, "\\gate_opt.csv"), gate_fid, ',')

