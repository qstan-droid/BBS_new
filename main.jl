include("functions//circuit.jl")
include("parameters.jl")

# want to vary two types of variables
# vary the order and block size
# vary the amplitude or the bias

# initialise the plots
samples_1 = []
samples_2 = []
fid_list = []
ave_fid = zeros(Float64, length(x))

for i = 1:length(x)
    if x_var == "alpha"
        if alpha_2 == false
            ave_fid[i], fid_list_temp, samples_1_temp, samples_2_temp = circuit(code, N_ord, [x[i], x[i]], block_size, err_place, err_info, measure, decode_type, sample_no, bias)
        else
            ave_fid[i], fid_list_temp, samples_1_temp, samples_2_temp = circuit(code, N_ord, [x[i], alpha_2], block_size, err_place, err_info, measure, decode_type, sample_no, bias)
        end
    elseif x_var == "bias"
        if where_bias == 1
            ave_fid[i], fid_list_temp, samples_1_temp, samples_2_temp = circuit(code, N_ord, alpha, block_size, err_place, err_info, measure, decode_type, sample_no, [x[i], 0])
        elseif where_bias == 2
            ave_fid[i], fid_list_temp, samples_1_temp, samples_2_temp = circuit(code, N_ord, alpha, block_size, err_place, err_info, measure, decode_type, sample_no, [0, x[i]])
        end
    end
    append!(samples_1, samples_1_temp)
    append!(samples_2, samples_2_temp)
    append!(fid_list, fid_list_temp)
end

# now we save onto a folder
