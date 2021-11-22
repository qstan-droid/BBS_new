# number of samples
sample_no = 100

# type of measurement
# heterodyne, opt_phase, adapt_homo
measure = ["opt_phase", "opt_phase"]

# note the invariable parameters of code
code = ["binomial", "binomial"]
block_size = [1, 3, 1] # row, col, rep
N_ord = [3, 3]

# define what errors are, how strong and where they're applied
# error_placement = [loss on block_1, dephase on block_1, loss on block_2, dephase on block_2]
# error_info = [nu_loss_1, nu_dephase_1, nu_loss_2, nu_dephase_2]
err_place = [true, false, true, false]
err_info = [0.1, 0.0, 0.1, 0.0]

# choose to vary alpha or bias
x_var = "alpha"

if x_var == "alpha"
    bias = [0, 0]

    x_min = 10
    x_step = 1
    x_max = 15

    x = x_min:x_step:x_max
    # same alpha for both blocks
    dif_alpha = false
    alpha_2 = 15
elseif x_var == "bias"
    alpha = [5, 15]
    where_bias = 2

    x_min = 0
    x_step = 0.1
    x_max = 2*pi

    x = x_min:x_step:x_max
end

# how to decode
# decoder_type = naive, naive_ave_bias, ml, ml_ave, mlnc
decode_type = "ml_ave"
err_spread_type = "normal" # normal, no_spread, err_trans