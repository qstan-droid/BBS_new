using QuantumOptics
using Distributions
include("measurement.jl")

function decoding(samples_1, samples_2, N_ord, block_size, err_place, err_info, sample_no, decode_type, measure, bias, xbasis)
    # initialise empty arrays
    outcomes_1 = zeros(Bool, sample_no)
    outcomes_2 = zeros(Bool, sample_no)

    # how to decode?
    if decode_type == "naive"

        # decode each qubit individually
        samples_out_1 = naive_decode(samples_1, N_ord[1], block_size, measure[1], 0)
        samples_out_2 = naive_decode(samples_2, N_ord[2], block_size, measure[2], 0)

        # decide outcome through these
        outcomes_1 = block_decode(samples_out_1, 1, block_size)
        outcomes_2 = block_decode(samples_out_2, 2, block_size)
    elseif decode_type == "bias"

        # decode each qubit individually
        samples_out_1 = naive_decode(samples_1, N_ord[1], block_size, measure[1], bias[1])
        samples_out_2 = naive_decode(samples_2, N_ord[2], block_size, measure[2], bias[2])

        # decide outcome through these
        outcomes_1 = block_decode(samples_out_1, 1, block_size)
        outcomes_2 = block_decode(samples_out_2, 2, block_size)
    elseif decode_type == "ave_max_like"

        # take the two samples and some information about the system and
        outcomes_1 = ave_max_like_decoder(samples_1)
        outcomes_2 = ave_max_like_decoder(samples_2)
    elseif decode_type == "max_likelihood"

        # Take the two samples and find the maximum likelihood
        outcomes_1, outcomes_2 = max_like_decoder(samples_1, samples_2, N_ord, err_info, xbasis, measure)
    end

    return outcomes_1, outcomes_2
end

#########################################################
# different functions, different decoding types

function naive_decode(samples, N_ord, block_size, measure, bias)
    # decide the outcome for each qubit individually
    row, col, sample_no = size(samples)
    samples_out = zeros(Int64, (row, col, sample_no))

    for k = 1:sample_no
        for i = 1:row
            for j = 1:col
                samples_out[i, j, k] = meas_outcome(samples[i, j, k], N_ord, measure, bias)
            end
        end
    end

    return samples_out
end

function max_like_decoder(samples_1, samples_2, N_ord, err_info, xbasis, measure)
    part = zeros(4)
    row, col, sample_no = size(samples_1)

    println(sample_no)

    xbasis_1 = xbasis[1]
    xbasis_2 = xbasis[2]

    N_ord_1 = N_ord[1]
    N_ord_2 = N_ord[2]

    # prepare measurement operators
    meas_op_1 = measurement_operator(measure[1], xbasis_1, N_ord_1)
    meas_op_2 = measurement_operator(measure[2], xbasis_2, N_ord_2)
    meas_ops = [meas_op_1, meas_op_2]

    # outcomes for majority
    maj_1 = zeros(Bool, (row, col, sample_no))
    maj_2 = zeros(Bool, (row, col, sample_no))

    # Fock space operators
    n_b_1 = xbasis_1[3]
    a_b_1 = xbasis_1[4]

    n_b_2 = xbasis_2[3]
    a_b_2 = xbasis_2[4]

    plus_1 = tensor(xbasis_1[1], dagger(xbasis_1[1]))
    min_1 = tensor(xbasis_1[2], dagger(xbasis_1[2]))

    plus_2 = tensor(xbasis_2[1], dagger(xbasis_2[1]))
    min_2 = tensor(xbasis_2[2], dagger(xbasis_2[2]))

    # unpack error information
    # err_info = nu_loss_1, nu_dephase_1, nu_loss_2, nu_dephase_2
    nu_loss_1 = err_info[1]
    nu_dephase_1 = err_info[2]
    nu_loss_2 = err_info[3]
    nu_dephase_2 = err_info[4]

    # Loss kraus operators
    if nu_loss_2 != 0
        A = function(p1, p2) 
            (((1-exp(-nu_loss_1))^(p1/2))/sqrt(factorial(big(p1))))*exp(-nu_loss_1*dense(n_b_1)/2)*a_b_1^p1 * exp(-1im*(p2*pi/(N_ord_1*N_ord_2))*dense(n_b_2))
        end
    else
        A = function(p1, p2) 
            (((1-exp(-nu_loss_1))^(p1/2))/sqrt(factorial(big(p1))))*exp(-nu_loss_1*dense(n_b_1)/2)*a_b_1^p1
        end
    end

    if nu_loss_1 != 0
        B = function(p1, p2) 
            #(((1-exp(-nu_loss_2))^(p1/2))/sqrt(factorial(big(p1))))*exp(-nu_loss_2*dense(n_b_2)/2)*a_b_2^p1 * exp(-1im*(p2*pi/(N_ord_1*N_ord_2))*dense(n_b_2))
            exp(-1im*(p2*pi/(N_ord_1*N_ord_2))*dense(n_b_2))
        end
    else
        B = function(p1, p2)
            (((1-exp(-nu_loss_2))^(p1/2))/sqrt(factorial(big(p1))))*exp(-nu_loss_2*dense(n_b_2)/2)*a_b_2^p1
        end
    end

    # find all likelihoods and choose the max

    for i = collect(1:sample_no)
        meas_ops_1 = meas_ops[1](samples_1[i])
        meas_ops_2 = meas_ops[2](samples_2[i])

        if nu_loss_1 == 0.0 && nu_loss_2 == 0.0
            part[1] = norm(tr(meas_ops_1*plus_1)*tr(meas_ops_2*plus_2))
            part[2] = norm(tr(meas_ops_1*plus_1)*tr(meas_ops_2*min_2))
            part[3] = norm(tr(meas_ops_1*min_1)*tr(meas_ops_2*plus_2))
            part[4] = norm(tr(meas_ops_1*min_1)*tr(meas_ops_2*min_2))
        else
            # precompute the traces
            A_plus = zeros(Complex{Float64}, 101)
            A_min = zeros(Complex{Float64}, 101)
            B_plus = zeros(Complex{Float64}, 101)
            B_min = zeros(Complex{Float64}, 101)

            for j = 0:100
                A_plus[j+1] = tr(meas_ops_1*(A(j, 0)*plus_1*dagger(A(j, 0))))
                A_min[j+1] = tr(meas_ops_1*(A(j, 0)*min_1*dagger(A(j, 0))))
                B_plus[j+1] = tr(meas_ops_2*(B(0, j)*plus_2*dagger(B(0, j))))
                B_min[j+1] = tr(meas_ops_2*(B(0, j)*min_2*dagger(B(0, j))))
            end

            part[1] = norm(sum(A_plus[j]*B_plus[j] for j = 1:101))
            part[2] = norm(sum(A_plus[j]*B_min[j] for j = 1:101))
            part[3] = norm(sum(A_min[j]*B_plus[j] for j = 1:101))
            part[4] = norm(sum(A_min[j]*B_min[j] for j = 1:101))
        end

        max_index = findmax(part)[2]

        if max_index == 1
            maj_1[i] = true
            maj_2[i] = true
        elseif max_index == 2
            maj_1[i] = true
            maj_2[i] = false
        elseif max_index == 3
            maj_1[i] = false
            maj_2[i] = true
        elseif max_index == 4
            maj_1[i] = false
            maj_2[i] = false
        end
    end

    return maj_1, maj_2
end

#########################################################
# functions which are tools for decoding

function meas_outcome(meas, N_ord, meas_type, bias)
    if meas_type == "heterodyne"
        phi = mod2pi(mod2pi(angle(meas)) + pi/(N_ord*2)) + bias
        if phi == 2*pi
            phi = 0
        end
        if phi > 2*pi
            phi = phi - 2*pi
        end
    elseif meas_type == "opt_phase"
        phi = mod2pi(mod2pi(convert(Float64, meas)) + pi/(N_ord*2)) + bias
        if phi == 2*pi
            phi = 0
        end
        if phi > 2*pi
            phi = phi - 2*pi
        end
    end

    k = 0
    edge = false

    while k < N_ord

        if phi > 0 + 2*k*pi/N_ord && phi < pi/N_ord + 2*k*pi/N_ord
            return +1
            edge = true
        elseif phi > pi/N_ord + 2*pi*k/N_ord && phi < (2*pi/N_ord) + 2*k*pi/N_ord
            return -1
            edge = true
        end
        k += 1
    end

    if edge == false
        coin = rand(1:2)
        if coin == 1
            return +1
        else
            return -1
        end
    end
end

function block_decode(samples_out, block_no, block_size)
    # unpack variables
    rep = block_size[3]
    row, col, sample_no = size(samples_out)
    maj = zeros(Bool, (row, col, sample_no))

    if block_no == 1
        # we compute total parity of each column
        # take majority vote over all parities
        for k = 1:sample_no
            col_par = fill(1, col)

            for i = 1:col
                col_par[i] = prod(samples_out[j, i, k] for j = 1:row)
            end

            # take majority vote over all column parities
            no_pos = count(l->(l == 1), col_par)
            no_neg = count(l->(l == -1), col_par)

            if no_pos > no_neg
                maj[k] = true
            else
                maj[k] = false
            end
        end

    elseif block_no == 2
        # take parities of each row
        for k = 1:sample_no
            row_par = fill(1, (rep, row, sample_no))

            # take row parities
            for k = 1:sample_no
                for i = 1:row
                    for j = 1:rep
                        row_par[j, i, k] = prod(samples_out[i, l, k] for l = (j-1)*div(col, rep)+1:j*div(col, rep))
                    end
                end
            end
            # majority vote over repetitions
            row_par_real = fill(1, (row, sample_no))
            for k = 1:sample_no
                for i = 1:row
                    par_maj = count(l->(l==1), row_par[:, i, k])
                    if par_maj > rep/2
                        row_par_real[i, k] = 1
                    else
                        row_par_real[i, k] = -1
                    end
                end
            end

            # majority vote over all column parities
            maj = zeros(Bool, sample_no)
            for k = 1:sample_no
                row_maj = count(l->(l==1), row_par_real[:, k])
                if row_maj > row/2
                    maj[k] = true
                else
                    maj[k] = false
                end
            end
        end
    end
    return maj
end
