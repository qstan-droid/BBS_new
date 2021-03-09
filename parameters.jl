# number of samples
sample_no = 1000

# type of measurement
measure = ["heterodyne", "heterodyne"]

# note the invariable parameters of code
code = ["binomial", "binomial"]
block_size = [1, 1, 1]
N_ord = [1, 1]

# define what errors are, how strong and where they're applied
err_place = [true, false, false, false]
err_info = [0.0, 0.0, 0.0, 0.0]

# choose to vary alpha or bias
x_var = "alpha"
if x_var == "alpha"
    bias = [0, 0]

    x_min = 2
    x_step = 2
    x_max = 14

    x = x_min:x_step:x_max
    # same alpha for both blocks
    dif_alpha = true
    alpha_2 = 50
elseif x_var == "bias"
    alpha = [10, 10]
    where_bias = 1

    x_min = 0
    x_step = 0.1
    x_max = 2*pi

    x = x_min:x_step:x_max
end

# how to decode
# decoder_type = naive, bias, ave_max_like, max_likelihood
decode_type = "naive"
