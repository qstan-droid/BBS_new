include("functions//decode.jl")

# make a random matrix (m by n matrix)
m = 3 # no of columns
n = 3 # no of rows
r = 3 # no of repetitions

# number of samples
k = 1

######## Testing block decode #########
#samples_output = rand([1, -1], (n, m*r, k))
block_size = [n, m, r]
#
#println("results input: ", samples_output)
#
## send to function to decode
#decoded_outcomes_1 = block_decode(samples_output, 1, block_size)
#decoded_outcomes_2 = block_decode(samples_output, 2, block_size)
#
#println("decoded output for block 1: ", decoded_outcomes_1)
#println("decoded output for block 2: ", decoded_outcomes_2)

######### Testing Error propagation #######
# Create a bunch of loss errors
#loss_1 = loss_sample(true, 0.1, block_size)