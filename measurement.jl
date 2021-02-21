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

# Finds expectation value of measurements
function meas_exp(meas_op, sample, err_prep_plus, err_prep_min)
    meas_exp_plus = dagger(err_prep_plus)*meas_op(sample)*err_prep_plus
    meas_exp_min = dagger(err_prep_min)*meas_op(sample)*err_prep_min
    meas_exp_pm = dagger(err_prep_plus)*meas_op(sample)*err_prep_min
    meas_exp_mp = dagger(err_prep_min)*meas_op(sample)*err_prep_plus

    return meas_exp_plus, meas_exp_min, meas_exp_pm, meas_exp_mp
end

###########################
function measurement_samples(err_prep_1, err_prep_2, err_exp_1, err_exp_2, block_no, measure_type, xbasis, N_ord, samples_1, code, block_size)

    #### Rejection sampling
    meas_op_1 = measurement_operator(measure_type[1], xbasis[1], N_ord[1])
    meas_op_2 = measurement_operator(measure_type[2], xbasis[2], N_ord[2])
    meas_ops = [meas_op_1, meas_op_2]

    # prepare the container for samples
    if block_no == 1
        row, col, sample_no = size(err_prep_1)
    elseif block_no == 2
        row, col, sample_no = size(err_prep_2)
    end

    samples = zeros(Complex(Float64), row, col, sample_no)
    norms = zeros(Complex(Float64), row, col, sample_no)

    meas_exp_plus = zeros(Complex(Float64), row, col, sample_no)
    meas_exp_min = zeros(Complex(Float64), row, col, sample_no)
    meas_exp_pm = zeros(Complex(Float64), row, col, sample_no)
    meas_exp_mp = zeros(Complex(Float64), row, col, sample_no)
    
    meas_exp = [meas_exp_plus, meas_exp_min, meas_exp_pm, meas_exp_mp]

    for k = 1:sample_no
        for i = 1:row
            for j = 1:col
                samples[i, j, k], meas_exp[1][i, j, k], meas_exp[2][i, j, k], meas_exp[3][i, j, k], meas_exp[4][i, j, k] = rejection_sampling(block_no, samples, samples_1, err_prep_1, err_prep_2, err_exp_1, err_exp_2, code, measure_type, [i, j, k], xbasis, meas_ops, meas_exp, block_size)
            end
        end
    end

    return samples, meas_exp[1], meas_exp[2], meas_exp[3], meas_exp[4]
end

function rejection_sampling(block_no, samples, samples_1, err_prep_1, err_prep_2, err_exp_1, err_exp_2, code, measure_type, loc, xbasis, meas_ops, meas_exp, block_size)

    if block_no == 1

        # find the envelope constant function
        ceil_constant = find_max_dist(block_no, samples , samples_1, measure_type, meas_ops, xbasis, code)*1.1
        counter = false

        while counter == false
            # sample a measurement
            samples[loc[1], loc[2], loc[3]] = sample_generator(code, measure_type, xbasis[1])
            meas_exp_plus[1][loc[1], loc[2], loc[3]], meas_exp_min[2][loc[1], loc[2], loc[3]], meas_exp_pm[2][loc[1], loc[2], loc[3]], meas_exp_mp[2][loc[1], loc[2], loc[3]] = meas_exp(samples[loc[1], loc[2], loc[3]], err_prep_1[1], err_prep_1[2])

            f_x = pdf_1(meas_exp[1][:, :, loc[3]], err_exp_1, err_exp_2, loc, block_size)

            # sample a random number between 1 and max
            u = rand(Uniform(0, ceil_constant))

            # check if condition is true
            if u < f_x
                counter = true
            end
        end

    elseif block_no == 2
        ceil_constant = find_max_dist(block_no)*1.1
        counter = false

        while counter == false
            # sample a measurement
            samples[loc[1], loc[2], loc[3]] = sample_generator(code, measure_type, xbasis[2])
            f_x = pdf_2()

            # sample a random number
            u = rand(Uniform(0, ceil_constant))

            # condition
            if u < f_x
                counter = true
            end
        end
    end

    norms[loc[1], loc[2], loc[3]] = f_x

    return samples[loc[1], loc[2], loc[3]], norms[loc[1], loc[2], loc[3]], meas_exp_plus[loc[1], loc[2], loc[3]], meas_exp_min[loc[1], loc[2], loc[3]]
end

function find_max_dist(block_no, samples, samples_1, meas_type, meas_ops, xbasis, code)

    #######
    if meas_type == "heterodyne"
        if code == "cat"
            edge = abs(xbasis[8])
        elseif code == "binomial"
            edge = abs(convert(Integer, round(sqrt(xbasis[8]))))
        end
        overflow = 0.5

        # create a range of values
        rad_range = 0:0.1:(edge + overflow)
        phi_range = 0:0.1:2*pi

        samples_range = zeros(length(rad_range), length(phi_range))

        for i = 1:length(rad_range)
            for j = 1:length(phi_range)
                samples_range[i, j] = rad_range[i, j]*(cos(phi_range[i, j]) + sin(phi_range[i, j])*1im)
            end
        end

    elseif meas_type == "opt_phase"
        samples_range = 0:0.1:2*pi
    end

    # find heights for every range and then choose the maximum
    if block_no == 1
        if meas_type == "heterodyne"
            row_range, col_range = size(samples_range)
            heights = zeros(row_range, col_range)

            for i = 1:row_range
                for j = 1:col_range
                    heights[i, j] = pdf_1()
                end
            end

        elseif meas_type == "opt_phase"
            heights = zeros(length(samples_range))

            for i = 1:length(samples_range)
                heights[i] = pdf_1()
            end
        end
    elseif block_no == 2
        if meas_type == "heterodyne"
            row_range, col_range = size(samples_range)
            heights = zeros(row_range, col_range)

            for i = 1:row_range
                for j = 1:col_range
                    heights[i, j] = pdf_2()
                end
            end

        elseif meas_type == "opt_phase"
            heights = zeros(length(samples_range))

            for i = 1:length(samples_range)
                heights[i] = pdf_2()
            end
        end
    end

    # return the maximum height
    max_height = findmax(heights)[1]
    return max_height
end

function sample_generator(code, meas_type, xbasis)

    if meas_type == "heterodyne"
        if code == "cat"
            edge = abs(xbasis[8])
        elseif code == "binomial"
            edge = abs(convert(Integer, round(sqrt(xbasis[8]))))
        end
        overflow = 3

        # sample over a circle for lesser surface area
        rad = rand(Uniform(0, edge + overflow))
        phi = rand(Uniform(0, 2*pi))
        beta = rad*cos(phi) + rad*sin(phi)*1im
    elseif meas_type == "opt_phase"
        beta = rand(Uniform(0, 2*pi))
    end

    return beta
end

#######################
# Probability functions
function pdf_1(meas_exp_1, err_exp_1, err_exp_2, loc, block_size)

    # setup all of the matrices
    no_row = block_size[1]
    no_col = block_size[2]
    rep = block_size[3]

    err_exp_2_zero = err_exp_2[1][:, :, loc[3]]
    err_exp_2_one = err_exp_2[2][:, :, loc[3]]
    err_exp_2_zo = err_exp_2[3][:, :, loc[3]]
    err_exp_2_oz = err_exp_2[4][:, :, loc[3]]

    err_exp_1_plus = err_exp_1[1][:, :, loc[3]]
    err_exp_1_min = err_exp_1[2][:, :, loc[3]]


    # pdf begins

    B = sum(sum( for j = 1:no_col) for i = 1:no_row)

end

function pdf_2()

end
