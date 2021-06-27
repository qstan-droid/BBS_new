# number of samples
sample_no = 1000

# type of measurement
measure = ["heterodyne", "heterodyne"]

# note the invariable parameters of code
code = ["binomial", "binomial"]
block_size = [1, 1, 1] # row, col, rep
N_ord = [3, 1]

# define what errors are, how strong and where they're applied
# error_placement = [loss on block_1, dephase on block_1, loss on block_2, dephase on block_2]
# error_info = [nu_loss_1, nu_dephase_1, nu_loss_2, nu_dephase_2]
err_place = [true, false, false, false]
err_info = [0.01, 0.0, 0.0, 0.0]

# choose to vary alpha or bias
x_var = "alpha"

if x_var == "alpha"
    bias = [0, 0]

    x_min = 1
    x_step = 1
    x_max = 15

    x = x_min:x_step:x_max
    # same alpha for both blocks
    dif_alpha = true
    alpha_2 = 15
elseif x_var == "bias"
    alpha = [10, 10]
    where_bias = 1

    x_min = 0
    x_step = 0.1
    x_max = 2*pi

    x = x_min:x_step:x_max
end

# how to decode
# decoder_type = naive, bias, ave_max_like, max_likelihood, max_likelihood_no_corr
decode_type = "max_likelihood"