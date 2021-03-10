using QuantumOptics
using Distributions

function measurement_operator(meas_type, xbasis, N_ord)

    coh_space = FockBasis(xbasis[9])

    if meas_type == "heterodyne"
        meas = function(x)
            tensor(coherentstate(coh_space, x), dagger(coherentstate(coh_space, x)))
        end

    elseif meas_type == "opt_phase"
        # Define the measurement operators
        meas = function(x)
            sum(exp(1im*n*x)*fockstate(coh_space, n) for n = 0:xbasis[8]*N_ord)/(xbasis[8]*N_ord + 1)
        end
    elseif meas_type == "homodyne"

    end

    return meas
end

# Finds expectation value of measurements
function meas_exp_prep(meas_op, sample, err_prep_plus, err_prep_min)

    meas_exp_plus = dagger(err_prep_plus)*meas_op(sample)*err_prep_plus
    meas_exp_min = dagger(err_prep_min)*meas_op(sample)*err_prep_min
    meas_exp_pm = dagger(err_prep_plus)*meas_op(sample)*err_prep_min
    meas_exp_mp = dagger(err_prep_min)*meas_op(sample)*err_prep_plus

    return meas_exp_plus, meas_exp_min, meas_exp_pm, meas_exp_mp
end

###########################
function measurement_samples(err_prep_1, err_prep_2, err_exp_1, err_exp_2, block_no, measure_type, meas_exp_1, xbasis, N_ord, samples_1, norms_1, code, block_size)

    #### Rejection sampling
    meas_op_1 = measurement_operator(measure_type[1], xbasis[1], N_ord[1])
    meas_op_2 = measurement_operator(measure_type[2], xbasis[2], N_ord[2])
    meas_ops = [meas_op_1, meas_op_2]

    # prepare the container for samples
    if block_no == 1
        row, col, sample_no = size(err_prep_1[1])
    elseif block_no == 2
        row, col, sample_no = size(err_prep_2[1])
    end

    samples = zeros(Complex{Float64}, row, col, sample_no)
    norms = fill(1.0+0.0*1im, row, col, sample_no)

    meas_exp_plus = zeros(Complex{Float64}, row, col, sample_no)
    meas_exp_min = zeros(Complex{Float64}, row, col, sample_no)
    meas_exp_pm = zeros(Complex{Float64}, row, col, sample_no)
    meas_exp_mp = zeros(Complex{Float64}, row, col, sample_no)

    meas_exp = [meas_exp_plus, meas_exp_min, meas_exp_pm, meas_exp_mp]

    for k = 1:sample_no
        for i = 1:row
            for j = 1:col
                samples[i, j, k], norms[i, j, k], meas_exp[1][i, j, k], meas_exp[2][i, j, k], meas_exp[3][i, j, k], meas_exp[4][i, j, k] = rejection_sampling(err_prep_1, err_prep_2, err_exp_1, err_exp_2, block_no, measure_type, xbasis, N_ord, samples, samples_1, code, block_size, meas_ops, meas_exp, [i, j, k], meas_exp_1, norms, norms_1)
            end
        end
    end

    return samples, norms, meas_exp[1], meas_exp[2], meas_exp[3], meas_exp[4]
end

function rejection_sampling(err_prep_1, err_prep_2, err_exp_1, err_exp_2, block_no, measure_type, xbasis, N_ord, samples, samples_1, code, block_size, meas_ops, meas_exp, loc, meas_exp_1, norms, norms_1)

    # find the envelope constant function
    # ceil_constant = find_max_dist(block_size, block_no, measure_type, meas_ops, err_prep_1, err_prep_2, err_exp_1, err_exp_2, meas_exp, meas_exp_1, code, loc)*1.1
    ceil_constant = 1
    counter = false

    # unpack loc
    x = loc[1]
    y = loc[2]
    z = loc[3]

    if block_no == 1
        while counter == false
            # sample a measurement
            samples[x, y, z] = sample_generator(code[1], measure_type[1], xbasis[1])
            meas_exp[1][x, y, z], meas_exp[2][x, y, z], meas_exp[3][x, y, z], meas_exp[4][x, y, z] = meas_exp_prep(meas_ops[block_no], samples[x, y, z], err_prep_1[1][x, y, z], err_prep_1[2][x, y, z])

            global f_x = pdf_1(meas_exp, err_exp_1, err_exp_2, norms, loc, block_size)

            # sample a random number between 1 and max
            u = rand(Uniform(0, ceil_constant))

            # check if condition is true
            println(abs(f_x))
            if u < abs(f_x)
                counter = true
            end
        end

    elseif block_no == 2
        while counter == false
            # sample a measurement
            samples[x, y, z] = sample_generator(code[2], measure_type[2], xbasis[2])
            meas_exp[1][x, y, z], meas_exp[2][x, y, z], meas_exp[3][x, y, z], meas_exp[4][x, y, z] = meas_exp_prep(meas_ops[block_no], samples[x, y, z], err_prep_2[1][x, y, z], err_prep_2[2][x, y, z])

            global f_x = pdf_2(meas_exp_1, meas_exp, err_exp_2, norms, norms_1, loc, block_size)

            # sample a random number
            u = rand(Uniform(0, ceil_constant))

            # condition
            if u < abs(f_x)
                counter = true
            end
        end
    end

    norms[x, y, z] = f_x

    return samples[x, y, z], norms[x, y, z], meas_exp[1][x, y, z], meas_exp[2][x, y, z], meas_exp[3][x, y, z], meas_exp[4][x, y, z]
end

function find_max_dist(block_size, block_no, meas_type, meas_ops, err_prep_1, err_prep_2, err_exp_1, err_exp_2, meas_exp, meas_exp_1, xbasis, code, loc)

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

        samples_range = zeros(Complex{Float64}, (length(rad_range), length(phi_range)))

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
                    meas_exp_1[1][loc[1], loc[2], loc[3]], meas_exp_1[2][loc[1], loc[2], loc[3]], meas_exp_1[3][loc[1], loc[2], loc[3]], meas_exp_1[4][loc[1], loc[2], loc[3]] = meas_exp(meas_ops[1], samples_range[i, j], err_prep_1[1], err_prep_1[2])
                    heights[i, j] = pdf_1(meas_exp_1, err_exp_1, err_exp_2, norms, loc, block_size)
                end
            end

        elseif meas_type == "opt_phase"
            heights = zeros(length(samples_range))

            for i = 1:length(samples_range)
                meas_exp_1[1][loc[1], loc[2], loc[3]], meas_exp_1[2][loc[1], loc[2], loc[3]], meas_exp_1[3][loc[1], loc[2], loc[3]], meas_exp_1[4][loc[1], loc[2], loc[3]] = meas_exp(meas_ops[1], samples_range[i, j], err_prep_1[1], err_prep_1[2])
                heights[i] = pdf_1(meas_exp_1, err_exp_1, err_exp_2, norms, loc, block_size)
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
function pdf_1(meas_exp_1, err_exp_1, err_exp_2, norms, loc, block_size)

    # setup all of the matrices
    no_row = block_size[1]
    no_col = block_size[2]
    rep = block_size[3]

    err_exp_1_plus = err_exp_1[1][:, :, loc[3]]
    err_exp_1_min = err_exp_1[2][:, :, loc[3]]
    err_exp_1_pm = err_exp_1[3][:, :, loc[3]]
    err_exp_1_mp = err_exp_1[4][:, :, loc[3]]

    err_exp_2_zero = err_exp_2[1][:, :, loc[3]]
    err_exp_2_one = err_exp_2[2][:, :, loc[3]]
    err_exp_2_zo = err_exp_2[3][:, :, loc[3]]
    err_exp_2_oz = err_exp_2[4][:, :, loc[3]]

    meas_exp_1_plus = meas_exp_1[1][:, :, loc[3]]
    meas_exp_1_min = meas_exp_1[2][:, :, loc[3]]
    meas_exp_1_pm = meas_exp_1[3][:, :, loc[3]]
    meas_exp_1_mp = meas_exp_1[4][:, :, loc[3]]

    norms = norms[:, :, loc[3]]

    ### pdf begins ###
    # the constant factor
    B_plus = prod(prod(err_exp_2_zero[i, j] for j = 1:no_col) + prod(err_exp_2_one[i, j] for j = 1:no_col) + prod(err_exp_2_zo[i, j] for j = 1:no_col) + prod(err_exp_2_oz[i, j] for j = 1:no_col) for i = 1:no_row*rep)/(2^(no_row*rep))
    B_min = prod(prod(err_exp_2_zero[i, j] for j = 1:no_col) + prod(err_exp_2_one[i, j] for j = 1:no_col) - prod(err_exp_2_zo[i, j] for j = 1:no_col) - prod(err_exp_2_oz[i, j] for j = 1:no_col) for i = 1:no_row*rep)/(2^(no_row*rep))
    B = B_plus + B_min

    # the variable factor
    lom = [0 0 0; 1 0 1; 0 1 1; 1 1 0] # location of minuses
    A = zeros(Complex{Float64}, 4)

    for a = 1:4
        if no_row == 1
            Before = 1
            After = 1
            if block_size[2] == 1
                Current = meas_exp_1_plus[loc[1], loc[2]] +
                            (-1)^lom[a, 1]*meas_exp_1_pm[loc[1], loc[2]] +
                            (-1)^lom[a, 2]*meas_exp_1_mp[loc[1], loc[2]] +
                            (-1)^lom[a, 3]*meas_exp_1_min[loc[1], loc[2]]
            else
                if loc[2] == no_col
                    Current = (prod(meas_exp_1_plus[loc[1], j] for j = 1:no_col) +
                                (-1)^lom[a, 1]*prod(meas_exp_1_pm[loc[1], j] for j = 1:no_col) +
                                (-1)^lom[a, 2]*prod(meas_exp_1_mp[loc[1], j] for j = 1:no_col) +
                                (-1)^lom[a, 3]*prod(meas_exp_1_min[loc[1], j] for j = 1:no_col))
                else
                    Current = (prod(meas_exp_1_plus[loc[1], j] for j = 1:loc[2])*prod(err_exp_1_plus[loc[1], j] for j = (loc[2] + 1):no_col) +
                                (-1)^lom[a, 1]*prod(meas_exp_1_pm[loc[1], j] for j = 1:loc[2])*prod(err_exp_1_pm[loc[1], j] for j = (loc[2] + 1):no_col) +
                                (-1)^lom[a, 2]*prod(meas_exp_1_mp[loc[1], j] for j = 1:loc[2])*prod(err_exp_1_mp[loc[1], j] for j = (loc[2] + 1):no_col) +
                                (-1)^lom[a, 3]*prod(meas_exp_1_min[loc[1], j] for j = 1:loc[2])*prod(err_exp_1_min[loc[1], j] for j = (loc[2] + 1):no_col))
                end
            end
        elseif loc[1] == 1
            Before = 1
            if loc[2] == no_col
                Current = (prod(meas_exp_1_plus[loc[1], j] for j = 1:no_col) +
                            (-1)^lom[a, 1]*prod(meas_exp_1_pm[loc[1], j] for j = 1:no_col) +
                            (-1)^lom[a, 2]*prod(meas_exp_1_mp[loc[1], j] for j = 1:no_col)+
                            (-1)^lom[a, 3]*prod(meas_exp_1_min[loc[1], j] for j = 1:no_col))
            else
                Current = (prod(meas_exp_1_plus[loc[1], j] for j = 1:loc[2])*prod(err_exp_1_plus[loc[1], j] for j = (loc[2] + 1):no_col) +
                            (-1)^lom[a, 1]*prod(meas_exp_1_pm[loc[1], j] for j = 1:loc[2])*prod(err_exp_1_pm[loc[1], j] for j = (loc[2] + 1):no_col) +
                            (-1)^lom[a, 2]*prod(meas_exp_1_mp[loc[1], j] for j = 1:loc[2])*prod(err_exp_1_mp[loc[1], j] for j = (loc[2] + 1):no_col) +
                            (-1)^lom[a, 3]*prod(meas_exp_1_min[loc[1], j] for j = 1:loc[2])*prod(err_exp_1_min[loc[1], j] for j = (loc[2] + 1):no_col))
            end
            After = prod(prod(err_exp_1_plus[i, j] for j = 1:no_col) +
                    (-1)^lom[a, 1]*prod(err_exp_1_pm[i, j] for j = 1:no_col) +
                    (-1)^lom[a, 2]*prod(err_exp_1_mp[i, j] for j = 1:no_col) +
                    (-1)^lom[a, 3]*prod(err_exp_1_min[i, j] for j = 1:no_col) for i = (loc[1]+1):no_row)
        elseif loc[1] == no_row
            Before = prod(prod(meas_exp_1_plus[i, j] for j = 1:no_col) +
                        (-1)^lom[a, 1]*prod(meas_exp_1_pm[i, j] for j = 1:no_col) +
                        (-1)^lom[a, 2]*prod(meas_exp_1_mp[i, j] for j = 1:no_col) +
                        (-1)^lom[a, 3]*prod(meas_exp_1_min[i, j] for j = 1:no_col) for i = 1:(loc[1]-1))
            if loc[2] == no_col
                Current = (prod(meas_exp_1_plus[loc[1], j] for j = 1:no_col) +
                            (-1)^lom[a, 1]*prod(meas_exp_1_pm[loc[1], j] for j = 1:no_col) +
                            (-1)^lom[a, 2]*prod(meas_exp_1_mp[loc[1], j] for j = 1:no_col) +
                            (-1)^lom[a, 3]*prod(meas_exp_1_min[loc[1], j] for j = 1:no_col))
            else
                Current = (prod(meas_exp_1_plus[loc[1], j] for j = 1:loc[2])*prod(err_exp_1_plus[loc[1], j] for j = (loc[2]+1):no_col) +
                            (-1)^lom[a, 1]*prod(meas_exp_1_pm[loc[1], j] for j = 1:loc[2])*prod(err_exp_1_plus[loc[1], j] for j = (loc[2]+1):no_col) +
                            (-1)^lom[a, 2]*prod(meas_exp_1_mp[loc[1], j] for j = 1:loc[2])*prod(err_exp_1_plus[loc[1], j] for j = (loc[2]+1):no_col) +
                            (-1)^lom[a, 3]*prod(meas_exp_1_min[loc[1], j] for j = 1:loc[2])*prod(err_exp_1_plus[loc[1], j] for j = (loc[2]+1):no_col))
            end
            After = 1
        else
            Before = prod(prod(meas_exp_1_plus[i, j] for j = 1:no_col) +
                        (-1)^lom[a, 1]*prod(meas_exp_1_pm[i, j] for j = 1:no_col) +
                        (-1)^lom[a, 2]*prod(meas_exp_1_mp[i, j] for j = 1:no_col) +
                        (-1)^lom[a, 3]*prod(meas_exp_1_min[i, j] for j = 1:no_col) for i = 1:(loc[1]-1))
            if loc[2] == no_col
                Current = (prod(meas_exp_1_plus[loc[1], j] for j = 1:no_col) +
                            (-1)^lom[a, 1]*prod(meas_exp_1_pm[loc[1], j] for j = 1:no_col) +
                            (-1)^lom[a, 2]*prod(meas_exp_1_mp[loc[1], j] for j = 1:no_col) +
                            (-1)^lom[a, 3]*prod(meas_exp_1_min[loc[1], j] for j = 1:no_col))
            else
                Current = (prod(meas_exp_1_plus[loc[1], j] for j = 1:loc[2])*prod(err_exp_1_plus[loc[1], j] for j = (loc[2]+1):no_col) +
                            (-1)^lom[a, 1]*prod(meas_exp_1_pm[loc[1], j] for j = 1:loc[2])*prod(err_exp_1_plus[loc[1], j] for j = (loc[2]+1):no_col) +
                            (-1)^lom[a, 2]*prod(meas_exp_1_mp[loc[1], j] for j = 1:loc[2])*prod(err_exp_1_plus[loc[1], j] for j = (loc[2]+1):no_col) +
                            (-1)^lom[a, 3]*prod(meas_exp_1_min[loc[1], j] for j = 1:loc[2])*prod(err_exp_1_plus[loc[1], j] for j = (loc[2]+1):no_col))
            end
            After = prod(prod(err_exp_1_plus[i, j] for j = 1:no_col) +
                    (-1)^lom[a, 1]*prod(err_exp_1_pm[i, j] for j = 1:no_col) +
                    (-1)^lom[a, 2]*prod(err_exp_1_mp[i, j] for j = 1:no_col) +
                    (-1)^lom[a, 3]*prod(err_exp_1_min[i, j] for j = 1:no_col) for i = (loc[1]+1):no_row)
        end

        A[a] = Before*Current*After
    end

    # product of all norm
    current_norm = 1
    for i = 1:loc[1]
        for j = 1:loc[2]
            current_norm = current_norm*norms[i, j]
        end
    end

    ans = B*sum(A[i] for i = 1:4)/(2^(no_row+1) * 4)
    return ans
end

function pdf_2(meas_exp_1, meas_exp_2, err_exp_2, norms, norms_1, loc, block_size)

    # setup all of the matrices
    no_row = block_size[1]
    no_col = block_size[2]
    rep = block_size[3]

    err_exp_2_zero = err_exp_2[1][:, :, loc[3]]
    err_exp_2_one = err_exp_2[2][:, :, loc[3]]
    err_exp_2_zo = err_exp_2[3][:, :, loc[3]]
    err_exp_2_oz = err_exp_2[4][:, :, loc[3]]

    # measurement expectation values
    meas_exp_1_plus = meas_exp_1[1][:, :, loc[3]]
    meas_exp_1_min = meas_exp_1[2][:, :, loc[3]]
    meas_exp_1_pm = meas_exp_1[3][:, :, loc[3]]
    meas_exp_1_mp = meas_exp_1[4][:, :, loc[3]]

    meas_exp_2_zero = meas_exp_2[1][:, :, loc[3]]
    meas_exp_2_one = meas_exp_2[2][:, :, loc[3]]
    meas_exp_2_zo = meas_exp_2[3][:, :, loc[3]]
    meas_exp_2_oz = meas_exp_2[4][:, :, loc[3]]

    norms = norms[:, :, loc[3]]
    norms_1 = norms_1[:, :, loc[3]]

    ### pdf begins ###
    # the constant
    lom = [0 0 0; 1 0 1; 0 1 1; 1 1 0] # location of minuses
    A = zeros(Complex{Float64}, 4)

    for a = 1:4
        A[a] = prod(prod(meas_exp_1_plus[i, j] for j = 1:no_col) +
                (-1)^lom[a, 1]*prod(meas_exp_1_pm[i, j] for j = 1:no_col) +
                (-1)^lom[a, 2]*prod(meas_exp_1_mp[i, j] for j = 1:no_col) +
                (-1)^lom[a, 3]*prod(meas_exp_1_min[i, j] for j = 1:no_col) for i = 1:no_row)
    end

    # get product of all block 1 norms
    current_norm_1 = 1
    for i = 1:no_row
        for j = 1:no_col
            current_norm_1 = current_norm_1*norms_1[i, j]
        end
    end

    A = sum(A[i] for i = 1:4)/(2^(no_row + 1) * current_norm_1)

    # the variable factor
    B = zeros(Complex{Float64}, 2)
    lom2 = [0 0 0; 1 1 0]

    for b = 1:2
        if no_row*rep == 1
            Before = 1
            After = 1
            if no_col == 1
                Current = meas_exp_2_zero[loc[1], loc[2]] +
                            (-1)^lom2[b, 1]*meas_exp_2_zo[loc[1], loc[2]] +
                            (-1)^lom2[b, 2]*meas_exp_2_oz[loc[1], loc[2]] +
                            (-1)^lom2[b, 3]*meas_exp_2_one[loc[1], loc[2]]
            else
                if loc[2] == no_col
                    Current = prod(meas_exp_2_zero[loc[1], j] for j = 1:no_col) +
                                (-1)^lom2[b, 1]*prod(meas_exp_2_zo[loc[1], j] for j = 1:no_col) +
                                (-1)^lom2[b, 2]*prod(meas_exp_2_oz[loc[1], j] for j = 1:no_col) +
                                (-1)^lom2[b, 3]*prod(meas_exp_2_one[loc[1], j] for j = 1:no_col)
                else
                    Current = prod(meas_exp_2_zero[loc[1], j] for j = 1:loc[2])*prod(err_exp_2_zero[loc[1], j] for j = loc[2]+1:no_col) +
                                (-1)^lom2[b, 1]*prod(meas_exp_2_zo[loc[1], j] for j = 1:no_col)*prod(err_exp_2_zo[loc[1], j] for j = loc[2]+1:no_col) +
                                (-1)^lom2[b, 2]*prod(meas_exp_2_oz[loc[1], j] for j = 1:no_col)*prod(err_exp_2_oz[loc[1], j] for j = loc[2]+1:no_col) +
                                (-1)^lom2[b, 3]*prod(meas_exp_2_one[loc[1], j] for j = 1:no_col)*prod(err_exp_2_one[loc[1], j] for j = loc[2]+1:no_col)
                end
            end
        elseif loc[1] == 1
            Before = 1
            if loc[2] == no_col
                Current = (prod(meas_exp_2_zero[loc[1], j] for j = 1:no_col) +
                            (-1)^lom2[b, 1]*prod(meas_exp_2_zo[loc[1], j] for j = 1:no_col) +
                            (-1)^lom2[b, 2]*prod(meas_exp_2_oz[loc[1], j] for j = 1:no_col) +
                            (-1)^lom2[b, 3]*prod(meas_exp_2_one[loc[1], j] for j = 1:no_col))
            else
                Current = (prod(meas_exp_2_zero[loc[1], j] for j = 1:loc[2])*prod(err_exp_2_zero[loc[1], j] for j = loc[2]+1:no_col) +
                            (-1)^lom2[b, 1]*prod(meas_exp_2_zo[loc[1], j] for j = 1:loc[2])*prod(err_exp_2_zo[loc[1], j] for j = loc[2]+1:no_col) +
                            (-1)^lom2[b, 2]*prod(meas_exp_2_oz[loc[1], j] for j = 1:loc[2])*prod(err_exp_2_oz[loc[1], j] for j = loc[2]+1:no_col) +
                            (-1)^lom2[b, 3]*prod(meas_exp_2_one[loc[1], j] for j = 1:loc[2])*prod(err_exp_2_one[loc[1], j] for j = loc[2]+1:no_col))
            end
            After = prod(prod(err_exp_2_zero[i, j] for j = 1:no_col) +
                        (-1)^lom2[b, 1]*prod(err_exp_2_zo[i, j] for j = 1:no_col) +
                        (-1)^lom2[b, 2]*prod(err_exp_2_oz[i, j] for j = 1:no_col) +
                        (-1)^lom2[b, 3]*prod(err_exp_2_one[i, j] for j = 1:no_col) for i = (loc[1] + 1):no_row*rep)
        elseif loc[1] == no_row
            Before = prod(prod(meas_exp_2_zero[i, j] for j = 1:no_col) +
                        (-1)^lom2[b, 1]*prod(meas_exp_2_zo[i, j] for j = 1:no_col) +
                        (-1)^lom2[b, 2]*prod(meas_exp_2_oz[i, j] for j = 1:no_col) +
                        (-1)^lom2[b, 3]*prod(meas_exp_2_one[i, j] for j = 1:no_col) for i = 1:(loc[1]-1))
            if loc[2] == no_col
                Current = (prod(meas_exp_2_zero[loc[1], j] for j = 1:no_col) +
                            (-1)^lom2[b, 1]*prod(meas_exp_2_zo[loc[1], j] for j = 1:no_col) +
                            (-1)^lom2[b, 2]*prod(meas_exp_2_oz[loc[1], j] for j = 1:no_col) +
                            (-1)^lom2[b, 3]*prod(meas_exp_2_one[loc[1], j] for j = 1:no_col))
            else
                Current = prod(meas_exp_2_zero[loc[1], j] for j = 1:loc[2])*prod(err_exp_2_zero[loc[1], j] for j = (loc[2]+1):no_col) +
                            (-1)^lom2[b, 1]*prod(meas_exp_2_zo[loc[1], j] for j = 1:loc[2])*prod(err_exp_2_zo[loc[1], j] for j = (loc[2]+1):no_col) +
                            (-1)^lom2[b, 2]*prod(meas_exp_2_oz[loc[1], j] for j = 1:loc[2])*prod(err_exp_2_oz[loc[1], j] for j = (loc[2]+1):no_col) +
                            (-1)^lom2[b, 3]*prod(meas_exp_2_one[loc[1], j] for j = 1:loc[2])*prod(err_exp_2_one[loc[1], j] for j = (loc[2]+1):no_col)
            end


            After = 1
        else
            Before = prod(prod(meas_exp_2_zero[i, j] for j = 1:no_col) +
                        (-1)^lom2[b, 1]*prod(meas_exp_2_zo[i, j] for j = 1:no_col) +
                        (-1)^lom2[b, 2]*prod(meas_exp_2_oz[i, j] for j = 1:no_col) +
                        (-1)^lom2[b, 3]*prod(meas_exp_2_one[i, j] for j = 1:no_col) for i = 1:(loc[1]-1))

            if loc[2] == no_col
                Current = (prod(meas_exp_2_zero[loc[1], j] for j = 1:no_col) +
                            (-1)^lom2[b, 1]*prod(meas_exp_2_zo[loc[1], j] for j = 1:no_col) +
                            (-1)^lom2[b, 2]*prod(meas_exp_2_oz[loc[1], j] for j = 1:no_col) +
                            (-1)^lom2[b, 3]*prod(meas_exp_2_one[loc[1], j] for j = 1:no_col))
            else
                Current = (prod(meas_exp_2_zero[loc[1], j] for j = 1:loc[2])*prod(err_exp_2_zero[loc[1], j] for j = (loc[2]+1):no_col) +
                            (-1)^lom2[b, 1]*prod(meas_exp_2_zo[loc[1], j] for j = 1:loc[2])*prod(err_exp_2_zo[loc[1], j] for j = (loc[2]+1):no_col) +
                            (-1)^lom2[b, 2]*prod(meas_exp_2_oz[loc[1], j] for j = 1:loc[2])*prod(err_exp_2_oz[loc[1], j] for j = (loc[2]+1):no_col) +
                            (-1)^lom2[b, 3]*prod(meas_exp_2_one[loc[1], j] for j = 1:loc[2])*prod(err_exp_2_one[loc[1], j] for j = (loc[2]+1):no_col))
            end

            After = prod(prod(err_exp_2_zero[i, j] for j = 1:no_col) +
                    (-1)^lom2[b, 1]*prod(err_exp_2_zo[i, j] for j = 1:no_col) +
                    (-1)^lom2[b, 2]*prod(err_exp_2_oz[i, j] for j = 1:no_col) +
                    (-1)^lom2[b, 3]*prod(err_exp_2_one[i, j] for j = 1:no_col) for i = (loc[1]+1):no_row)

        end

        B[b] = Before*Current*After
    end

    # get product of all block 2 norms
    current_norm_2 = 1
    for i = 1:loc[1]
        for j = 1:loc[2]
            current_norm_2 = current_norm_2*norms[i, j]
        end
    end


    ans = A*sum(B[i] for i = 1:2)/(2^(no_row) * 4 * current_norm_2)
    return ans
end
