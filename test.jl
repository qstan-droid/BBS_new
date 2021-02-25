include("functions//circuit.jl")

# parameters
code = ["binomial", "binomial"]
N_ord = [3, 3]
alpha = [10, 10]
block_size = [1, 1, 1]
err_place = [true, false, true, false]
err_info = [0.01, 0.0, 0.01, 0.0]
measure = ["heterodyne", "heterodyne"] # only heterodyne
decode_type = "naive" # naive or bias
sample_no = 10
bias = [0, 0]

if code == "cat"
    dim_1 = convert(Int64, round(2*alpha[1]^2 + alpha[1], digits=0))
    dim_2 = convert(Int64, round(2*alpha[2]^2 + alpha[2], digits=0))
elseif code == "binomial"
    dim_1 = (alpha[1]+1)*(N_ord[1])
    dim_2 = (alpha[2]+1)*(N_ord[2])
end

dim = [dim_1, dim_2]


@time ave_fid, fid_list, samples_1, samples_2 = circuit(code, N_ord, dim, alpha, block_size, err_place, err_info, measure, decode_type, sample_no, bias)
