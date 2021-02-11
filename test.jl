include("errors.jl")

# place to test functions
a = fill(1.0, (2, 3))
a = [1 1 1; 2 2 2]
b = fill(1.0, (2, 9))
c = fill(0.0, (2, 3))
d = fill(0.0, (2, 9))

a1, c1, b1, d1 = error_propagation(a, c, b, d, [2, 3, 3], [1, 1])
println("block_1_loss: ", a1)
println("block_1_dephase: ", c1)
println("block_2_loss: ", b1)
println("block_2_dephase: ", d1)
