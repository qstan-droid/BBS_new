using DelimitedFiles
include("functions//circuit.jl")
include("test_functions//test_functions.jl")
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

    println("-------------------------------------")
    println("x_vary: ", x_var, " | ", "x: ", x[i], " | order: ", N_ord )

    if x_var == "alpha"
        if dif_alpha == false
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

#args = [1, "naive_1"]

# now we save onto a folder
open("parameters.txt", "w") do file
#open(string("data_", ARGS[2], "/parameters.txt"), "w") do file
    write(file, string("code: ", code, "\n"))
    write(file, string("sample_no: ", sample_no, "\n"))
    write(file, string("n: ", block_size[1], "\n"))
    write(file, string("m: ", block_size[2], "\n"))
    write(file, string("rep: ", block_size[3], "\n"))
    write(file, string("N_ord: ", N_ord, "\n"))
    write(file, "-------------------------------------\n")
    write(file, string("loss_1: ", err_place[1], ", ", err_info[1], "\n"))
    write(file, string("loss_2: ", err_place[3], ", ", err_info[3], "\n"))
    write(file, string("dephase_1: ", err_place[2], ", ", err_info[2], "\n"))
    write(file, string("dephase_2: ", err_place[4], ", ", err_info[4], "\n"))
    write(file, "-------------------------------------\n")
    write(file, string("measurement_type: ", measure, "\n"))
    write(file, string("decode_type: ", decode_type, "\n"))
    write(file, "-------------------------------------\n")
    write(file, string("x_var: ", x_var, "\n"))
    if x_var == "alpha"
        write(file, string("alpha: ", x_min, "-", x_step, "-", x_max, "\n"))
        write(file, string("bias: ", bias, "\n"))
        if dif_alpha
            write(file, string("differ_alpha: ", dif_alpha, ", ", alpha_2, "\n"))
        end
    elseif x_var == "bias"
        write(file, string("bias: ", x_min, "-", x_step, "-", x_max, "\n"))
        write(file, string("alpha: ", alpha, "\n"))
    end

end

<<<<<<< HEAD
spot = 1
=======


spot = 4
>>>>>>> 57b717bbe118e3e4a621c20701257942cf42be6e
p1_plus, p1_min = samples_plot(samples_1[sample_no*(spot - 1) + 1:sample_no*spot], N_ord[1], x[spot], measure[1], bias[1])
p2_plus, p2_min = samples_plot(samples_2[sample_no*(spot - 1) + 1:sample_no*spot], N_ord[2], alpha_2, measure[2], bias[2])

#writedlm(string("data_", ARGS[2], "/", ARGS[1], ".csv"), ave_fid, ',')
savefig(p1_plus, "samples_1_plot_plus")
savefig(p1_min, "samples_1_plot_min")
savefig(p2_plus, "samples_2_plot_plus")
savefig(p2_min, "samples_2_plot_min")
