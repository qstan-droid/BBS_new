

# given a set of complex numbers, return complex plane picture
function samples_plot(samples, N_ord)
    row, col, sample_no = size(samples)


    # split into x and y, choose a particular qubit in a block
    x = zeros(Float64, sample_no)
    y = zeros(Float64, sample_no)

    for k = 1:sample_no
        x[k] = real(samples[1, 1, k])
        y[k] = imag(samples[1, 1, k])
    end


end
