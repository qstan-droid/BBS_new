using QuantumOptics
using Distributions

function measurement_operator(meas_type, xbasis, N_ord, H_mn)

    coh_space = FockBasis(xbasis[9])

    # Define the measurement operators
    if meas_type == "heterodyne"
        meas = function(x)
            tensor(coherentstate(coh_space, x), dagger(coherentstate(coh_space, x)))/pi
        end

    elseif meas_type == "opt_phase"
        
        meas = function(x)
            # sum(exp(1im*n*x)*fockstate(coh_space, n) for n = 0:xbasis[8]*N_ord)/(xbasis[8]*N_ord + 1)
            sum(sum(exp(1im*(m-n)*x)*tensor(fockstate(coh_space, m), dagger(fockstate(coh_space, n))) for n = 0:xbasis[8]*N_ord) for m = 0:xbasis[8]*N_ord)/(xbasis[8]*N_ord + 1)
        end
    elseif meas_type == "adapt_homo"
        meas = function(x)
            sum(sum(exp(1im*x*(m-n))*H_mn[m+1, n+1]*tensor(fockstate(coh_space, m), dagger(fockstate(coh_space, n))) for n = 0:xbasis[8]*N_ord) for m = 0:xbasis[8]*N_ord)/(2*pi)
        end
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

#######################
function measurement_samples(err_prep_1, err_prep_2, err_exp_1, err_exp_2, block_no, measure_type, meas_exp_1, xbasis, N_ord, samples_1, norms_1, code, block_size, loss_norm_1, loss_norm_2, loss_1, loss_2, H_mn, alpha, err_info)

    #### Rejection sampling
    meas_op_1 = measurement_operator(measure_type[1], xbasis[1], N_ord[1], H_mn)
    meas_op_2 = measurement_operator(measure_type[2], xbasis[2], N_ord[2], H_mn)
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

    # acceptance fraction
    no_of_times_list = zeros(Complex{Float64}, row, col, sample_no)

    for k = 1:sample_no
        for i = 1:row
            for j = 1:col
                samples[i, j, k], norms[i, j, k], meas_exp[1][i, j, k], meas_exp[2][i, j, k], meas_exp[3][i, j, k], meas_exp[4][i, j, k], no_of_times_list[i, j, k] = rejection_sampling(err_prep_1, err_prep_2, err_exp_1, err_exp_2, block_no, measure_type, xbasis, N_ord, samples, samples_1, code, block_size, meas_ops, meas_exp, [i, j, k], meas_exp_1, norms, norms_1, loss_norm_1[:, :, k], loss_norm_2[:, :, k], loss_1[:, :, k], loss_2[:, :, k], alpha, err_info)
            end
        end
    end

    return samples, norms, meas_exp[1], meas_exp[2], meas_exp[3], meas_exp[4], no_of_times_list
end

function rejection_sampling(err_prep_1, err_prep_2, err_exp_1, err_exp_2, block_no, measure_type, xbasis, N_ord, samples, samples_1, code, block_size, meas_ops, meas_exp, loc, meas_exp_1, norms, norms_1, loss_norm_1, loss_norm_2, loss_1, loss_2, alpha, err_info)

    # find the envelope constant function
    # if it is single mode, use the constant ceiling
    # if it is the rep code, use non constant

    row = block_size[1]
    col = block_size[2]

    if row == 1 && col == 1
        ceil_constant = 0.6
    else
        #ceil_constant = abs(find_max_dist(block_size, block_no, measure_type[block_no], meas_ops, err_prep_1, err_prep_2, err_exp_1, err_exp_2, meas_exp, meas_exp_1, xbasis[block_no], code[block_no], loc, loss_norm_1, loss_norm_2, norms, norms_1))*1.1
        ceil_constant = 1.3
    end
    
    counter = false

    # unpack loc
    x = loc[1]
    y = loc[2]
    z = loc[3]

    no_of_times = 0
    #if 1 in loss
    #    println(loc)
    #    println(ceil_constant)
    #    println(loss)
    #    println(err_exp_1[1][:, :, z])
    #    println(meas_exp[1][:, :, z])
    #end

    if block_no == 1
        while counter == false
            # sample a measurement
            samples[x, y, z] = sample_generator(code[1], measure_type[1], xbasis[1])
            meas_exp[1][x, y, z], meas_exp[2][x, y, z], meas_exp[3][x, y, z], meas_exp[4][x, y, z] = meas_exp_prep(meas_ops[block_no], samples[x, y, z], err_prep_1[1][x, y, z], err_prep_1[2][x, y, z])

            global f_x = pdf_1(meas_exp, err_exp_1, err_exp_2, norms, loc, block_size, loss_norm_1, loss_norm_2)
            # global f_x = pdf_alt_1(code, N_ord, alpha, measure_type, err_spread_type, loss_1, loss_2, samples, norms, loss_norm_1, loss_norm_2, loc, err_info)

            if abs(f_x) == 0.0
                println("f_x is 0")
            end

            # sample a random number between 1 and max
            u = rand(Uniform(0, ceil_constant))

            # check if condition is true
            if abs(f_x) > ceil_constant
                println("goes above one: ", abs(f_x))
            else
                if u < abs(f_x)
                    no_of_times += 1
                    counter = true
                else 
                    no_of_times += 1
                end
            end
        end

    elseif block_no == 2
        while counter == false
            # sample a measurement
            samples[x, y, z] = sample_generator(code[2], measure_type[2], xbasis[2])
            meas_exp[1][x, y, z], meas_exp[2][x, y, z], meas_exp[3][x, y, z], meas_exp[4][x, y, z] = meas_exp_prep(meas_ops[block_no], samples[x, y, z], err_prep_2[1][x, y, z], err_prep_2[2][x, y, z])

            global f_x = pdf_2(meas_exp_1, meas_exp, err_exp_2, norms, norms_1, loc, block_size, loss_norm_1, loss_norm_2)

            # sample a random number
            u = rand(Uniform(0, ceil_constant))

            # condition
            if abs(f_x) > ceil_constant
                println("goes above one: ", abs(f_x))
            else
                if u < abs(f_x)
                    no_of_times += 1
                    counter = true
                else
                    #println(abs(f_x))
                    no_of_times += 1
                end
            end
        end
    end

    norms[x, y, z] = f_x

    return samples[x, y, z], norms[x, y, z], meas_exp[1][x, y, z], meas_exp[2][x, y, z], meas_exp[3][x, y, z], meas_exp[4][x, y, z], 1/no_of_times
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
    elseif meas_type == "opt_phase" || meas_type == "adapt_homo"
        beta = rand(Uniform(0, 2*pi))
    end

    return beta
end

#######################
# Probability functions
function pdf_1(meas_exp_1, err_exp_1, err_exp_2, norms, loc, block_size, loss_norm_1, loss_norm_2)

    # setup all of the matrices
    no_row = block_size[1]
    no_col = block_size[2]
    rep = block_size[3]

    err_exp_1_zero = err_exp_1[1][:, :, loc[3]]
    err_exp_1_one = err_exp_1[2][:, :, loc[3]]
    err_exp_1_zo = err_exp_1[3][:, :, loc[3]]
    err_exp_1_oz = err_exp_1[4][:, :, loc[3]]

    err_exp_2_zero = err_exp_2[1][:, :, loc[3]]
    err_exp_2_one = err_exp_2[2][:, :, loc[3]]
    err_exp_2_zo = err_exp_2[3][:, :, loc[3]]
    err_exp_2_oz = err_exp_2[4][:, :, loc[3]]

    meas_exp_1_zero = meas_exp_1[1][:, :, loc[3]]
    meas_exp_1_one = meas_exp_1[2][:, :, loc[3]]
    meas_exp_1_zo = meas_exp_1[3][:, :, loc[3]]
    meas_exp_1_oz = meas_exp_1[4][:, :, loc[3]]

    norms = norms[:, :, loc[3]]
    x = loc[1] # corresponds to row_number
    y = loc[2] # corresponds to col number
    z = loc[3]

    ### pdf begins ###
    ### B_factor ###
    B_plus = prod(prod(prod(err_exp_2_zero[i, (p-1)*no_col + j] for j = 1:no_col) + prod(err_exp_2_zo[i, (p-1)*no_col + j] for j = 1:no_col) + prod(err_exp_2_oz[i, (p-1)*no_col + j] for j = 1:no_col) + prod(err_exp_2_one[i, (p-1)*no_col + j] for j = 1:no_col) for p = 1:rep) for i = 1:no_row)/(2^(rep*no_row))
    B_min = prod(prod(prod(err_exp_2_zero[i, (p-1)*no_col + j] for j = 1:no_col) - prod(err_exp_2_zo[i, (p-1)*no_col + j] for j = 1:no_col) - prod(err_exp_2_oz[i, (p-1)*no_col + j] for j = 1:no_col) + prod(err_exp_2_one[i, (p-1)*no_col + j] for j = 1:no_col) for p = 1:rep) for i = 1:no_row)/(2^(rep*no_row))

    B_fac = B_plus + B_min

    ### A_factor ###
    #A_plus = prod(prod(meas_exp_1_zero[i, j] for i = 1:x)*prod(err_exp_1_zero[i, j] for i = x+1:no_row; init=1) + prod(meas_exp_1_zo[i, j] for i = 1:x)*prod(err_exp_1_zo[i, j] for i = x+1:no_row; init=1) + prod(meas_exp_1_oz[i, j] for i = 1:x)*prod(err_exp_1_oz[i, j] for i = x+1:no_row; init=1) + prod(meas_exp_1_one[i, j] for i = 1:x)*prod(err_exp_1_one[i, j] for i = x+1:no_row; init=1) for j = 1:y)*
    #            prod(prod(err_exp_1_zero[i, j] for i = 1:no_row) + prod(err_exp_1_zo[i, j] for i = 1:no_row) + prod(err_exp_1_oz[i, j] for i = 1:no_row) + prod(err_exp_1_one[i, j] for i = 1:no_row) for j = y+1:no_col; init = 1)/(2^no_col)
    #A_min = prod(prod(meas_exp_1_zero[i, j] for i = 1:x)*prod(err_exp_1_zero[i, j] for i = x+1:no_row; init=1) - prod(meas_exp_1_zo[i, j] for i = 1:x)*prod(err_exp_1_zo[i, j] for i = x+1:no_row; init=1) - prod(meas_exp_1_oz[i, j] for i = 1:x)*prod(err_exp_1_oz[i, j] for i = x+1:no_row; init=1) + prod(meas_exp_1_one[i, j] for i = 1:x)*prod(err_exp_1_one[i, j] for i = x+1:no_row; init=1) for j = 1:y)*
    #            prod(prod(err_exp_1_zero[i, j] for i = 1:no_row) - prod(err_exp_1_zo[i, j] for i = 1:no_row) - prod(err_exp_1_oz[i, j] for i = 1:no_row) + prod(err_exp_1_one[i, j] for i = 1:no_row) for j = y+1:no_col; init=1)/(2^no_col)
    
    A_plus = prod(prod(meas_exp_1_zero[i, j] for i = 1:x; init=1)*prod(err_exp_1_zero[i, j] for i = x+1:no_row; init=1) + prod(meas_exp_1_zo[i, j] for i = 1:x; init=1)*prod(err_exp_1_zo[i, j] for i = x+1:no_row; init=1) + prod(meas_exp_1_oz[i, j] for i = 1:x; init=1)*prod(err_exp_1_oz[i, j] for i = x+1:no_row; init=1) + prod(meas_exp_1_one[i, j] for i = 1:x; init=1)*prod(err_exp_1_one[i, j] for i = x+1:no_row; init=1) for j = 1:y; init=1)*
                prod(prod(meas_exp_1_zero[i, j] for i = 1:x-1; init=1)*prod(err_exp_1_zero[i, j] for i = x:no_row; init=1) + prod(meas_exp_1_zo[i, j] for i = 1:x-1; init=1)*prod(err_exp_1_zo[i, j] for i = x:no_row; init=1) + prod(meas_exp_1_oz[i, j] for i = 1:x-1; init=1)*prod(err_exp_1_oz[i, j] for i = x:no_row; init=1) + prod(meas_exp_1_one[i, j] for i = 1:x-1; init=1)*prod(err_exp_1_one[i, j] for i = x:no_row; init=1) for j = y+1:no_col; init=1)/(2^no_col)
    A_min = prod(prod(meas_exp_1_zero[i, j] for i = 1:x)*prod(err_exp_1_zero[i, j] for i = x+1:no_row; init=1) - prod(meas_exp_1_zo[i, j] for i = 1:x)*prod(err_exp_1_zo[i, j] for i = x+1:no_row; init=1) - prod(meas_exp_1_oz[i, j] for i = 1:x)*prod(err_exp_1_oz[i, j] for i = x+1:no_row; init=1) + prod(meas_exp_1_one[i, j] for i = 1:x)*prod(err_exp_1_one[i, j] for i = x+1:no_row; init=1) for j = 1:y)*
                prod(prod(meas_exp_1_zero[i, j] for i = 1:x-1; init=1)*prod(err_exp_1_zero[i, j] for i = x:no_row; init=1) - prod(meas_exp_1_zo[i, j] for i = 1:x-1; init=1)*prod(err_exp_1_zo[i, j] for i = x:no_row; init=1) - prod(meas_exp_1_oz[i, j] for i = 1:x-1; init=1)*prod(err_exp_1_oz[i, j] for i = x:no_row; init=1) + prod(meas_exp_1_one[i, j] for i = 1:x-1; init=1)*prod(err_exp_1_one[i, j] for i = x:no_row; init=1) for j = y+1:no_col; init=1)/(2^no_col)

    A_fac = A_plus + A_min

    ### norms ###
    # get product of all block 1 norms
    current_norm = 1
    for i = 1:no_row
        for j = 1:no_col
            current_norm = current_norm*norms[i, j]*loss_norm_1[i, j]*loss_norm_2[i, j]
        end
    end

    ### answer ###
    ans = A_fac*B_fac/(4*current_norm)

    return ans
end

function pdf_2(meas_exp_1, meas_exp_2, err_exp_2, norms, norms_1, loc, block_size, loss_norm_1, loss_norm_2)

    # setup all of the matrices
    no_row = block_size[1]
    no_col = block_size[2]
    rep = block_size[3]

    # location of current qubit
    x = loc[1] # row number
    y = loc[2] # col number
    z = loc[3]

    rho = Int(ceil(y/no_col)) # repetition number
    y_tru = Int(y - (rho - 1)*no_col) # row number in the repetition

    err_exp_2_zero = err_exp_2[1][:, :, loc[3]]
    err_exp_2_one = err_exp_2[2][:, :, loc[3]]
    err_exp_2_zo = err_exp_2[3][:, :, loc[3]]
    err_exp_2_oz = err_exp_2[4][:, :, loc[3]]

    # measurement expectation values
    meas_exp_1_zero = meas_exp_1[1][:, :, loc[3]]
    meas_exp_1_one = meas_exp_1[2][:, :, loc[3]]
    meas_exp_1_zo = meas_exp_1[3][:, :, loc[3]]
    meas_exp_1_oz = meas_exp_1[4][:, :, loc[3]]

    meas_exp_2_zero = meas_exp_2[1][:, :, loc[3]]
    meas_exp_2_one = meas_exp_2[2][:, :, loc[3]]
    meas_exp_2_zo = meas_exp_2[3][:, :, loc[3]]
    meas_exp_2_oz = meas_exp_2[4][:, :, loc[3]]

    norms = norms[:, :, loc[3]]
    norms_1 = norms_1[:, :, loc[3]]

    ### pdf begins ###
    ### A_factor ###
    A_plus = prod(prod(meas_exp_1_zero[i, j] for i = 1:no_row) + prod(meas_exp_1_zo[i, j] for i = 1:no_row) + prod(meas_exp_1_oz[i, j] for i = 1:no_row) + prod(meas_exp_1_one[i, j] for i = 1:no_row) for j = 1:no_col)/(2^no_col)
    A_min = prod(prod(meas_exp_1_zero[i, j] for i = 1:no_row) - prod(meas_exp_1_zo[i, j] for i = 1:no_row) - prod(meas_exp_1_oz[i, j] for i = 1:no_row) + prod(meas_exp_1_one[i, j] for i = 1:no_row) for j = 1:no_col)/(2^no_col)

    A_fac = A_plus + A_min

    ### B_factor ###
    B_plus = prod(prod(prod(meas_exp_2_zero[i, (p-1)*no_col + j] for j = 1:no_col; init=1) + prod(meas_exp_2_zo[i, (p-1)*no_col + j] for j = 1:no_col; init=1) + prod(meas_exp_2_oz[i, (p-1)*no_col + j] for j = 1:no_col; init=1) + prod(meas_exp_2_one[i, (p-1)*no_col + j] for j = 1:no_col; init=1) for p = 1:rep; init=1) for i = 1:x-1; init=1)*
                prod(prod(meas_exp_2_zero[x, (p-1)*no_col + j] for j = 1:no_col; init=1) + prod(meas_exp_2_zo[x, (p-1)*no_col + j] for j = 1:no_col; init=1) + prod(meas_exp_2_oz[x, (p-1)*no_col + j] for j = 1:no_col; init=1) + prod(meas_exp_2_one[x, (p-1)*no_col + j] for j = 1:no_col; init=1) for p = 1:rho-1; init=1)*
                (prod(meas_exp_2_zero[x, (rho-1)*no_col + j] for j = 1:y_tru; init=1)*prod(err_exp_2_zero[x, (rho-1)*no_col + j] for j = y_tru + 1:no_col; init=1) + prod(meas_exp_2_zo[x, (rho-1)*no_col + j] for j = 1:y_tru; init=1)*prod(err_exp_2_zo[x, (rho-1)*no_col + j] for j = y_tru + 1:no_col; init=1) + prod(meas_exp_2_oz[x, (rho-1)*no_col + j] for j = 1:y_tru; init=1)*prod(err_exp_2_oz[x, (rho-1)*no_col + j] for j = y_tru + 1:no_col; init=1) + prod(meas_exp_2_one[x, (rho-1)*no_col + j] for j = 1:y_tru; init=1)*prod(err_exp_2_one[x, (rho-1)*no_col + j] for j = y_tru+1:no_col; init=1))*
                prod(prod(err_exp_2_zero[x, (p-1)*no_col+j] for j = 1:no_col; init=1) + prod(err_exp_2_zo[x, (p-1)*no_col + j] for j = 1:no_col; init=1) + prod(err_exp_2_oz[x, (p-1)*no_col + j] for j = 1:no_col; init=1) + prod(err_exp_2_one[x, (p-1)*no_col + j] for j = 1:no_col; init=1) for p = rho + 1:rep; init=1)*
                prod(prod(prod(err_exp_2_zero[i, (p-1)*no_col+j] for j = 1:no_col; init=1) + prod(err_exp_2_zo[i, (p-1)*no_col+j] for j = 1:no_col; init=1) + prod(err_exp_2_oz[i, (p-1)*no_col+j] for j = 1:no_col; init=1) + prod(err_exp_2_one[i, (p-1)*no_col+j] for j = 1:no_col; init=1) for p = 1:rep; init=1) for i = x+1:no_row; init=1)/(2^(no_row*rep))
    
    B_min = prod(prod(prod(meas_exp_2_zero[i, (p-1)*no_col + j] for j = 1:no_col; init=1) - prod(meas_exp_2_zo[i, (p-1)*no_col + j] for j = 1:no_col; init=1) - prod(meas_exp_2_oz[i, (p-1)*no_col + j] for j = 1:no_col; init=1) + prod(meas_exp_2_one[i, (p-1)*no_col + j] for j = 1:no_col; init=1) for p = 1:rep; init=1) for i = 1:x-1; init=1)*
                prod(prod(meas_exp_2_zero[x, (p-1)*no_col + j] for j = 1:no_col; init=1) - prod(meas_exp_2_zo[x, (p-1)*no_col + j] for j = 1:no_col; init=1) - prod(meas_exp_2_oz[x, (p-1)*no_col + j] for j = 1:no_col; init=1) + prod(meas_exp_2_one[x, (p-1)*no_col + j] for j = 1:no_col; init=1) for p = 1:rho-1; init=1)*
                (prod(meas_exp_2_zero[x, (rho-1)*no_col + j] for j = 1:y_tru; init=1)*prod(err_exp_2_zero[x, (rho-1)*no_col + j] for j = y_tru + 1:no_col; init=1) - prod(meas_exp_2_zo[x, (rho-1)*no_col + j] for j = 1:y_tru; init=1)*prod(err_exp_2_zo[x, (rho-1)*no_col + j] for j = y_tru + 1:no_col; init=1) - prod(meas_exp_2_oz[x, (rho-1)*no_col + j] for j = 1:y_tru; init=1)*prod(err_exp_2_oz[x, (rho-1)*no_col + j] for j = y_tru + 1:no_col; init=1) + prod(meas_exp_2_one[x, (rho-1)*no_col + j] for j = 1:y_tru; init=1)*prod(err_exp_2_one[x, (rho-1)*no_col + j] for j = y_tru+1:no_col; init=1))*
                prod(prod(err_exp_2_zero[x, (p-1)*no_col+j] for j = 1:no_col; init=1) - prod(err_exp_2_zo[x, (p-1)*no_col + j] for j = 1:no_col; init=1) - prod(err_exp_2_oz[x, (p-1)*no_col + j] for j = 1:no_col; init=1) + prod(err_exp_2_one[x, (p-1)*no_col + j] for j = 1:no_col; init=1) for p = rho + 1:rep; init=1)*
                prod(prod(prod(err_exp_2_zero[i, (p-1)*no_col+j] for j = 1:no_col; init=1) - prod(err_exp_2_zo[i, (p-1)*no_col+j] for j = 1:no_col; init=1) - prod(err_exp_2_oz[i, (p-1)*no_col+j] for j = 1:no_col; init=1) + prod(err_exp_2_one[i, (p-1)*no_col+j] for j = 1:no_col; init=1) for p = 1:rep; init=1) for i = x+1:no_row; init=1)/(2^(no_row*rep))

    B_fac = B_plus + B_min

    ### norms ###
    current_norm = 1
    for i = 1:no_row
        for j = 1:no_col
            current_norm = current_norm*norms_1[i, j]*norms[i, j]*loss_norm_1[i, j]*loss_norm_2[i, j]
        end
    end

    ### answer ###
    ans = A_fac*B_fac/(4*current_norm)
    return ans
end

function pdf_alt_1(code, N_ord, alpha, measure_type, err_spread_type, loss_1, loss_2, samples, norms, loss_norm_1, loss_norm_2, loc, err_info)
    # prepare the beginning
    no_row = block_size[1]
    no_col = block_size[2]
    rep = block_size[3]

    # prepare code parameters
    code_1 = code[1] # the code of each block
    code_2 = code[2]

    N_ord_1 = N_ord[1]
    N_ord_2 = N_ord[2]

    alpha_1 = alpha[1]
    alpha_2 = alpha[2]

    meas_1 = measure_type[1] # type of measurement
    meas_2 = measure_type[2]

    # position and info
    x = loc[1] # row
    y = loc[2] # col
    z = loc[3]

    nu_1 = err_info[1]
    nu_2 = err_info[3]

    # prepare the loss and samples
    l_1 = loss_1
    l_2 = loss_2

    phi = samples[:, :, loc[3]]

    # noirms
    norms = norms[:, :, loc[3]]

    #---------------------------------
    # functions that we'll need
    G = function(k::Int64, K::Int64)
        sqrt(binomial(k, K)/(2^K))
    end

    c_nu = function(nu, l::Int64)
        ((1 - exp(-nu))^(l/2))/sqrt(factorial(big(l)))
    end

    P = function(n, r)
        factorial(big(n))/factorial(big(n - r))
    end

    del = function(n1, n2)
        if n1 == n2
            return 1
        else
            return 0
        end
    end
    #---------------------------------
    # begin writing the pdfs

    if err_spread_type == "no_spread"

        # for block 1
        A = function(K, N1, N2, nu, l, l2, w, p)
            (c_nu(nu, l)^2/2)*sum(sum((1 + (-1)^w *(-1)^k1)*(1 + (-1)^p *(-1)^k2)*G(k1, K)*G(k2, K)*exp(-nu*((k1 + k2)*N1 - 2*l)/2)*P(k1*N1, l)*P(k2*N1, l)*del(k1*N1 - l, k2*N1 - l) for k2=0:K) for k1=0:K)
        end

        if code_1 == "binomial"
            if meas_1 == "heterodyne"
                AM = function(K, N1, N2, nu, l, l2, x_samp, w, p)
                    (c_nu(nu, l)^2 /(2*pi))*sum(sum((1 + (-1)^w *(-1)^k1)*(1 + (-1)^p *(-1)^k2)*G(k1, K)*G(k2, K)*exp(-nu*((k1 + k2)*N1 - 2*l)/2)*P(k1*N1, l)*P(k2*N1, l)*exp(-abs(x_samp)^2)*(conj(x_samp)^(k1*N1 - l) /sqrt(factorial(big(k1*N1 - l))))*(x_samp^(k1*N1 - l) /sqrt(factorial(big(k2*N1 - l)))) for k1=0:K) for k2=0:K) 
                end
            elseif meas_1 == "opt_phase"
                AM = function(K, N1, N2, nu, l, l2, x_samp, w, p)
                    (c_nu(nu, l)^2 /(2*(K*N1+1)))*sum(sum((1 + (-1)^w *(-1)^k1)*(1 + (-1)^p *(-1)^k2)*G(k1, K)*G(k2, K)*exp(-nu*((k1 + k2)*N1 - 2*l)/2)*P(k1*N1, l)*P(k2*N1, l)*exp(1im*(k2-k1)*N1*x_samp) for k1=0:K) for k2=0:K) 
                end
            end
        end

        # for block 2
        if code_2 == "binomial"
            B = function(K, N1, N2, nu, l, l2, w, p)
                (c_nu(nu, l)^2/2)*sum(sum((1 + (-1)^w *(-1)^k1)*(1 + (-1)^p *(-1)^k2)*G(k1, K)*G(k2, K)*exp(-nu*((k1 + k2)*N2 - 2*l)/2)*P(k1*N2, l)*P(k2*N2, l)*del(k1*N2 - l, k2*N2 - l) for k2=0:K) for k1=0:K)
            end
        end

    elseif err_spread_type == "normal"

        # for block 1
        A = function(K, N1, N2, nu, l1, l2, w, p)
            (c_nu(nu, l1)/2)*sum(sum((1 + (-1)^w *(-1)^k1)*(1 + (-1)^p *(-1)^k2)*G(k1, K)*G(k2, K)*exp(-(k1*N1 - l1)*(nu/2 + 1im*l2*pi/(N1*N2)))*exp(-(k2*N1 - l1)*(nu/2 - 1im*l2*pi/(N1*N2)))*P(k1*N1, l1)*P(k2*N1, l1)*del(k1*N1 - l1, k2*N1 - l1) for k1=0:K) for k2=0:K)
        end

        if code_1 == "binomial"
            if meas_1 == "heterodyne"
                AM = function(K, N1, N2, nu, l1, l2, x_samp, w, p)
                    (c_nu(nu, l1)/(2*pi))*sum(sum((1 + (-1)^w *(-1)^k1)*(1 + (-1)^p *(-1)^k2)*G(k1, K)*G(k2, K)*exp(-(k1*N1 - l1)*(nu/2 + 1im*l2*pi/(N1*N2)))*exp(-(k2*N1 - l1)*(nu/2 - 1im*l2*pi/(N1*N2)))*P(k1*N1, l1)*P(k2*N1, l1)*exp(-abs(x_samp)^2)*(conj(x_samp)^(k1*N1 - l) /sqrt(factorial(big(k1*N1 - l))))*(x_samp^(k1*N1 - l) /sqrt(factorial(big(k2*N1 - l)))) for k1=0:K) for k2=0:K)
                end
            elseif meas_1 == "opt_phase"
                AM = function(K, N1, N2, nu, l1, l2, x_samp, w, p)
                    (c_nu(nu, l)^2 /(2*(K*N_ord_1+1)))*sum(sum((1 + (-1)^w *(-1)^k1)*(1 + (-1)^p *(-1)^k2)*G(k1, K)*G(k2, K)*exp(-(k1*N1 - l1)*(nu/2 + 1im*l2*pi/(N1*N2)))*exp(-(k2*N1 - l1)*(nu/2 - 1im*l2*pi/(N1*N2)))*P(k1*N1, l1)*P(k2*N1, l1)*exp(1im*(k2-k1)*N1*x_samp) for k1=0:K) for k2=0:K)
                end
            end
        end
        
        # for block 2
        if code_2 == "binomial"
            B = function(K, N1, N2, nu, l1, l2, w, p)
                (c_nu(nu, l2)/2)*sum(sum((1 + (-1)^w *(-1)^k1)*(1 + (-1)^p *(-1)^k2)*G(k1, K)*G(k2, K)*exp(-1im*l1*pi*k1*N2/(N1*N2))*exp(1im*l1*pi*k2*N2/(N1*N2))*exp(-nu*((k1+k2)*N2 - 2*l2)/2)*P(k1*N2, l2)*P(k2*N2, l2)*del(k1*N2 - l2, k2*N2 - l2) for k1=0:K) for k2=0:K)
            end
        end
    end

    # putting it all together

    # putting B together
    B_plus = prod(prod((prod(B(alpha_2, N_ord_1, N_ord_2, nu_2, l_1[i, j], l_2[i, (p-1)*no_col + j], 0, 0) for j = 1:no_col) + prod(B(alpha_1, N_ord_1, N_ord_2, nu_2, l_1[i, j], l_2[i, (p-1)*no_col + j], 0, 1) for j = 1:no_col) + prod(B(alpha_1, N_ord_1, N_ord_2, nu_2, l_1[i, j], l_2[i, (p-1)*no_col + j], 1, 0) for j = 1:no_col) + prod(B(alpha_1, N_ord_1, N_ord_2, nu_2, l_1[i, j], l_2[i, (p-1)*no_col + j], 1, 1) for j = 1:no_col)) for p = 1:rep) for i = 1:no_row)/(2^(no_row*rep))
    B_min = prod(prod((prod(B(alpha_2, N_ord_1, N_ord_2, nu_2, l_1[i, j], l_2[i, (p-1)*no_col + j], 0, 0) for j = 1:no_col) - prod(B(alpha_1, N_ord_1, N_ord_2, nu_2, l_1[i, j], l_2[i, (p-1)*no_col + j], 0, 1) for j = 1:no_col) - prod(B(alpha_1, N_ord_1, N_ord_2, nu_2, l_1[i, j], l_2[i, (p-1)*no_col + j], 1, 0) for j = 1:no_col) + prod(B(alpha_1, N_ord_1, N_ord_2, nu_2, l_1[i, j], l_2[i, (p-1)*no_col + j], 1, 1) for j = 1:no_col)) for p = 1:rep) for i = 1:no_row)/(2^(no_row*rep))
    B_factor = B_plus + B_min

    # putting A together
    A_plus = prod(prod(AM(alpha_1, N_ord_1, N_ord_2, nu_1, l_1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), phi[i, j], 0, 0)  for i=1:x)*prod(A(alpha_1, N_ord_1, N_ord_2, nu_1, l1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), 0, 0) for i=x+1:no_row;init=1) + prod(AM(alpha_1, N_ord_1, N_ord_2, nu_1, l_1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), phi[i, j], 0, 1)  for i=1:x)*prod(A(alpha_1, N_ord_1, N_ord_2, nu_1, l1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), 0, 1) for i=x+1:no_row;init=1) + prod(AM(alpha_1, N_ord_1, N_ord_2, nu_1, l_1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), phi[i, j], 1, 0)  for i=1:x)*prod(A(alpha_1, N_ord_1, N_ord_2, nu_1, l1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), 1, 0) for i=x+1:no_row;init=1) + prod(AM(alpha_1, N_ord_1, N_ord_2, nu_1, l_1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), phi[i, j], 1, 1)  for i=1:x)*prod(A(alpha_1, N_ord_1, N_ord_2, nu_1, l1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), 1, 1) for i=x+1:no_row;init=1) for j = 1:y)*
            prod(prod(AM(alpha_1, N_ord_1, N_ord_2, nu_1, l_1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), phi[i, j], 0, 0)  for i=1:x-1;init=1)*prod(A(alpha_1, N_ord_1, N_ord_2, nu_1, l1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), 0, 0) for i=x:no_row) + prod(AM(alpha_1, N_ord_1, N_ord_2, nu_1, l_1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), phi[i, j], 0, 1)  for i=1:x-1;init=1)*prod(A(alpha_1, N_ord_1, N_ord_2, nu_1, l1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), 0, 1) for i=x:no_row) + prod(AM(alpha_1, N_ord_1, N_ord_2, nu_1, l_1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), phi[i, j], 1, 0)  for i=1:x-1;init=1)*prod(A(alpha_1, N_ord_1, N_ord_2, nu_1, l1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), 1, 0) for i=x:no_row) + prod(AM(alpha_1, N_ord_1, N_ord_2, nu_1, l_1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), phi[i, j], 1, 1)  for i=1:x-1;init=1)*prod(A(alpha_1, N_ord_1, N_ord_2, nu_1, l1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), 1, 1) for i=x:no_row) for j=y+1:no_col;init=1)
    A_min = prod(prod(AM(alpha_1, N_ord_1, N_ord_2, nu_1, l_1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), phi[i, j], 0, 0)  for i=1:x)*prod(A(alpha_1, N_ord_1, N_ord_2, nu_1, l1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), 0, 0) for i=x+1:no_row;init=1) - prod(AM(alpha_1, N_ord_1, N_ord_2, nu_1, l_1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), phi[i, j], 0, 1)  for i=1:x)*prod(A(alpha_1, N_ord_1, N_ord_2, nu_1, l1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), 0, 1) for i=x+1:no_row;init=1) - prod(AM(alpha_1, N_ord_1, N_ord_2, nu_1, l_1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), phi[i, j], 1, 0)  for i=1:x)*prod(A(alpha_1, N_ord_1, N_ord_2, nu_1, l1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), 1, 0) for i=x+1:no_row;init=1) + prod(AM(alpha_1, N_ord_1, N_ord_2, nu_1, l_1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), phi[i, j], 1, 1)  for i=1:x)*prod(A(alpha_1, N_ord_1, N_ord_2, nu_1, l1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), 1, 1) for i=x+1:no_row;init=1) for j = 1:y)*
            prod(prod(AM(alpha_1, N_ord_1, N_ord_2, nu_1, l_1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), phi[i, j], 0, 0)  for i=1:x-1;init=1)*prod(A(alpha_1, N_ord_1, N_ord_2, nu_1, l1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), 0, 0) for i=x:no_row) - prod(AM(alpha_1, N_ord_1, N_ord_2, nu_1, l_1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), phi[i, j], 0, 1)  for i=1:x-1;init=1)*prod(A(alpha_1, N_ord_1, N_ord_2, nu_1, l1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), 0, 1) for i=x:no_row) - prod(AM(alpha_1, N_ord_1, N_ord_2, nu_1, l_1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), phi[i, j], 1, 0)  for i=1:x-1;init=1)*prod(A(alpha_1, N_ord_1, N_ord_2, nu_1, l1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), 1, 0) for i=x:no_row) + prod(AM(alpha_1, N_ord_1, N_ord_2, nu_1, l_1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), phi[i, j], 1, 1)  for i=1:x-1;init=1)*prod(A(alpha_1, N_ord_1, N_ord_2, nu_1, l1[i, j], sum(l_2[i, (p-1)*no_col + j] for p=1:rep), 1, 1) for i=x:no_row) for j=y+1:no_col;init=1)
    A_factor = A_plus + A_min

    # get the norms
    current_norm = 1
    for i = 1:no_row
        for j = 1:no_col
            current_norm = current_norm*norms[i, j]*loss_norm_1[i, j]*loss_norm_2[i, j]
        end
    end

    ## answer ##
    ans = A_factor*B_factor/(4*current_norm)

    return ans
        
end

#######################

function find_max_dist(block_size, block_no, meas_type, meas_ops, err_prep_1, err_prep_2, err_exp_1, err_exp_2, meas_exp, meas_exp_1, xbasis, code, loc, loss_norm_1, loss_norm_2, norms, norms_1)

    #######
    # prepare a samples range
    if meas_type == "heterodyne"
        if code == "cat"
            edge = abs(xbasis[8])
        elseif code == "binomial"
            edge = abs(convert(Integer, round(sqrt(xbasis[8]))))
        end
        overflow = 1.0

        # create a range of values
        rad_range = 0:0.1:(edge + overflow)
        phi_range = 0:0.1:2*pi

        samples_range = zeros(Complex{Float64}, (length(rad_range), length(phi_range)))

        for i = 1:length(rad_range)
            for j = 1:length(phi_range)
                samples_range[i, j] = rad_range[i]*(cos(phi_range[j]) + sin(phi_range[j])*1im)
            end
        end

    elseif meas_type == "opt_phase"
        samples_range = 0:0.1:2*pi
    end

    # find heights for every range and then choose the maximum
    if block_no == 1
        if meas_type == "heterodyne"
            row_range, col_range = size(samples_range)
            global heights = zeros(Complex{Float64}, (row_range, col_range))

            for i = 1:row_range
                for j = 1:col_range
                    meas_exp[1][loc[1], loc[2], loc[3]], meas_exp[2][loc[1], loc[2], loc[3]], meas_exp[3][loc[1], loc[2], loc[3]], meas_exp[4][loc[1], loc[2], loc[3]] = meas_exp_prep(meas_ops[1], samples_range[i, j], err_prep_1[1][loc[1], loc[2], loc[3]], err_prep_1[2][loc[1], loc[2], loc[3]])
                    heights[i, j] = pdf_1(meas_exp, err_exp_1, err_exp_2, norms, loc, block_size, loss_norm_1, loss_norm_2)
                end
            end
        elseif meas_type == "opt_phase"
            global heights = zeros(Complex{Float64}, length(samples_range))

            for i = 1:length(samples_range)
                meas_exp[1][loc[1], loc[2], loc[3]], meas_exp[2][loc[1], loc[2], loc[3]], meas_exp[3][loc[1], loc[2], loc[3]], meas_exp[4][loc[1], loc[2], loc[3]] = meas_exp_prep(meas_ops[1], samples_range[i], err_prep_1[1][loc[1], loc[2], loc[3]], err_prep_1[2][loc[1], loc[2], loc[3]])
                heights[i] = pdf_1(meas_exp, err_exp_1, err_exp_2, norms, loc, block_size, loss_norm_1, loss_norm_2)
            end
        end
    elseif block_no == 2
        if meas_type == "heterodyne"
            row_range, col_range = size(samples_range)
            global heights = zeros(Complex{Float64}, (row_range, col_range))

            for i = 1:row_range
                for j = 1:col_range
                    meas_exp[1][loc[1], loc[2], loc[3]], meas_exp[2][loc[1], loc[2], loc[3]], meas_exp[3][loc[1], loc[2], loc[3]], meas_exp[4][loc[1], loc[2], loc[3]] = meas_exp_prep(meas_ops[2], samples_range[i, j], err_prep_2[1][loc[1], loc[2], loc[3]], err_prep_2[2][loc[1], loc[2], loc[3]])
                    heights[i, j] = pdf_2(meas_exp_1, meas_exp, err_exp_2, norms, norms_1, loc, block_size, loss_norm_1, loss_norm_2)
                end
            end

        elseif meas_type == "opt_phase"
            global heights = zeros(Complex{Float64}, length(samples_range))

            for i = 1:length(samples_range)
                meas_exp[1][loc[1], loc[2], loc[3]], meas_exp[2][loc[1], loc[2], loc[3]], meas_exp[3][loc[1], loc[2], loc[3]], meas_exp[4][loc[1], loc[2], loc[3]] = meas_exp_prep(meas_ops[2], samples_range[i], err_prep_2[1][loc[1], loc[2], loc[3]], err_prep_2[2][loc[1], loc[2], loc[3]])
                heights[i] = pdf_2(meas_exp_1, meas_exp, err_exp_2, norms, norms_1, loc, block_size, loss_norm_1, loss_norm_2)
            end
        end
    end
    
    # return the maximum height
    max_height = findmax(abs.(heights))[1]
    return max_height
end