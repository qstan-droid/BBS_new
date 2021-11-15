using Base: find_source_file
using DelimitedFiles

# which folder
parent_folder = "rep mode stuff"
folder = string(parent_folder, "\\data_1_3_1_ml_ave_0.01_same_err_no_diff_N2_opt_more_samp")

# collecting all the files
no_of_files = 10
samples_no = 500

# the means
fid_mean = readdlm(string(@__DIR__, "\\", folder, "\\opt.csv"), ',', Float64)
gate_fid_mean = readdlm(string(@__DIR__, "\\", folder, "\\gate_opt.csv"), ',', Float64)

# Calculate the SD
fid_SD = zeros(length(fid_mean))
gate_fid_SD = zeros(length(gate_fid_mean))

#for i = 1:length(fid_mean)-1
for i = 1:11
    for j = 1:no_of_files
        fid_list = readdlm(string(@__DIR__, "\\", folder, "\\fidelity_list_", j, ".csv"), ',', Float64)
        gate_fid_list = readdlm(string(@__DIR__, "\\", folder, "\\gate_fidelity_list_", j, ".csv"), ',', Float64)

        fid_SD[i] = fid_SD[i] + sum((fid_list[i, k] - fid_mean[i])^2 for k = 1:samples_no)
        gate_fid_SD[i] = gate_fid_SD[i] + sum((gate_fid_list[i, k] - gate_fid_mean[i])^2 for k = 1:samples_no)
    end
end


fid_SD = sqrt.(fid_SD ./ samples_no*no_of_files)
gate_fid_SD = sqrt.(gate_fid_SD ./ samples_no*no_of_files)

fid_SE = fid_SD ./ samples_no
gate_fid_SE = gate_fid_SD ./ samples_no

writedlm(string(@__DIR__, "\\", folder, "\\SE_opt.csv"), fid_SE, ',')
writedlm(string(@__DIR__, "\\", folder, "\\gate_SE_opt.csv"), gate_fid_SE, ',')