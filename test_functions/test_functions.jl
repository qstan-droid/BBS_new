using Plots

# given a set of complex numbers, return complex plane picture
function samples_plot(samples, N_ord, alpha, meas_type, bias)

    # split into x and y, choose a particular qubit in a block
    x_plus = zeros(Float64, sample_no)
    y_plus = zeros(Float64, sample_no)

    x_min = zeros(Float64, sample_no)
    y_min = zeros(Float64, sample_no)

    for k = 1:sample_no
        verdict = meas_outcome(samples[k], N_ord, meas_type, bias)

        if verdict == 1
            x_plus[k] = real(samples[k])
            y_plus[k] = imag(samples[k])
        elseif verdict == -1
            x_min[k] = real(samples[k])
            y_min[k] = imag(samples[k])
        end

    end



    p_plus = scatter(x_plus, y_plus, xlim = (-alpha-3, alpha+3), ylim = (-alpha-3, alpha+3))
    p_min = scatter(x_min, y_min, xlim = (-alpha-3, alpha+3), ylim = (-alpha-3, alpha+3))

    return p_plus, p_min
end
