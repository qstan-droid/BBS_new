using QuantumOptics
using Distributions

function measuremet_operator(meas_type, xbasis, N_ord)

    if meas_type == "heterodyne"
        meas = function(x)
            coherentstate(xbasis[9], x)/sqrt(pi)
        end

    elseif meas_type == "opt_phase"
        # Define the measurement operators
        meas = function(x)
            sum(exp(1im*n*x)*fockstate(xbasis[9], n) for n = 0:xbasis[8]*N_ord)/(xbasis[8]*N_ord + 1)
        end
    elseif meas_type == "homodyne"

    end

    return meas
end

###########################
function measurement_samples(err_prep_1, err_prep_2, block_no, measure_type, xbasis, N_ord)

    #### Rejection sampling
    meas_op_1 = measurement_operator(measure_type[1], xbasis[1], N_ord[1])
    meas_op_2 = measurement_operator(measure_type[2], xbasis[2], N_ord[2])

    row, col, sample_no = size(err_prep_1)

end
