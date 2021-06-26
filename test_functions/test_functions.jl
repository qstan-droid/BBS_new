using Plots
using StatsPlots

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

function plot_samples(samples, norms, measure, block_no)
    no_row, no_col, sample_size = size(samples)

    # split samples into x and y
    x = zeros(no_row, no_col, sample_size)
    y = zeros(no_row, no_col, sample_size)
    norm_abs = zeros(no_row, no_col, sample_size)
    for i = 1:sample_size
        for r = 1:no_row
            for c = 1:no_col
                if measure == "heterodyne"
                    x[r, c, i] = real(samples[r, c, i])
                    y[r, c, i] = imag(samples[r, c, i])
                elseif measure == "opt_phase"
                    x[r, c, i] = cos(real(samples[r, c, i]))
                    y[r, c, i] = sin(real(samples[r, c, i]))
                end
                norm_abs[r, c, i] = abs(norms[r, c, i])
            end
        end
    end

    # plot shit
    for r = 1:no_row
        for c = 1:no_col
            plot3D = scatter(x[r, c, :], y[r, c, :], norm_abs[r, c, :])
            plotPlane = scatter(x[r, c, :], y[r, c, :])

            savefig(plot3D, string("plot3D_", measure, "_(", r, ",", c, ")_", block_no))
            savefig(plotPlane, string("plotPlane_", measure, "_(", r, ",", c, ")_", block_no))
        end
    end
end

####################################################
# TEST PLOTTING #

function initial_plotting(ave_fid, ave_fid_SE, x)
    
    p = plot(x, 1 .- ave_fid, yerr = ave_fid_SE, yaxis = :log)
    display(p)
    return p
end

function loss_hist(loss)
    max = findmax(loss)[1]
    hist = zeros(Int64, max+1)

    for i = 0:max
        hist[i+1] = count(j->(j == i), loss)
    end
    display(histogram(vec(loss)))
    println(vec(hist) ./ 1000)
end