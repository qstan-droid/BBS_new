# number of samples
sample_no = 1

# type of measurement
measure = ["heterodyne", "heterodyne"]

# note the invariable parameters of code
code = ["binomial", "binomial"]
block_size = [3, 3, 1]
N_ord = [1, 1]

# define what errors are, how strong and where they're applied
err_place = [true, false, true, false]
err_info = [0.0, 0.0, 0.0, 0.0]

# choose to vary alpha or bias
x_var = "alpha"
if x_var == "alpha"
    bias = [0, 0]
    x = 2:2:2
    # same alpha for both blocks
    dif_alpha = true
    alpha_2 = 50
elseif x_var == "bias"
    alpha = [10, 10]
    where_bias = 1
    x = 0:0.1:2*pi
end

# how to decode
decode_type = "naive"
