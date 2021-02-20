using Distributions

function loss_sample(err_place, nu_loss, m, n, xbasis, sample_no)
    loss_samples = zeros(n, m, sample_no)
    loss_norm = fill(1.0, (n, m, sample_no))

    if err_place == true
        for k = 1:sample_no
            for i = 1:n
                for j = 1:m
                    loss_samples[i, j, k], loss_norm[i, j, k] = loss_sample(nu_loss, xbasis)
                end
            end
        end
    end

    return loss_samples, loss_norm
end

function dephase_sample(err_place, nu_dephase, m, n, sample_no)
    dephase_norm = fill(1.0, (n, m, sample_no))
    d = Normal(0, nu_dephase)

    dephase_samples = zeros(n, m, sample_no)

    if err_place == true
        dephase_samples = rand(d, (n, m, sample_no))

        for k = 1:sample_no
            for i = 1:n
                for j = 1:m
                    dephase_norm[i, j, k] = exp((dephase_samples[i, j, k]/nu_dephase)^2 /2)/(nu_dephase*sqrt(2*pi))
                end
            end
        end
    end

    return dephase_samples, dephase_norm
end

function error_propagation(loss_1, dephase_1, loss_2, dephase_2, block, N_ord)
    n = block[1]
    m = block[2]
    no_rep = block[3]

    for i = 1:n
        # do the spread from block 1 to block 2
        for j = 1:m
            for k = 1:no_rep
                dephase_2[i, m*(k - 1)+j] = dephase_2[i, m*(k - 1)+j] + pi*loss_1[i, j]/(N_ord[1]*N_ord[2])
            end
        end

        # then do spread from block 2 to block 1
        for j = 1:m
            for k = 1:no_rep
                dephase_1[i, j] = dephase_1[i, j] + pi*loss_2[i, m*(k - 1)+j]/(N_ord[1]*N_ord[2])
            end
        end
    end

    return loss_1, dephase_1, loss_2, dephase_2
end

#### ERROR SAMPLING ####

function loss_pdf(nu_l, k_l, xbasis)
    A = (((1 - exp(-nu_l))^(k_l/2))/sqrt(factorial(big(k_l))))*exp(-nu_l*dense(xbasis[3])/2)*xbasis[4]^k_l

    ans = (dagger(A*xbasis[1])*A*xbasis[1] + dagger(A*xbasis[2])*A*xbasis[2])/2
    return ans
end

function discrete_loss_sampling(nu_l, xbasis)
    # initialise
    k::Int64 = 0
    s = loss_pdf(nu_l, k, xbasis)
    n = rand(Unifor(0, 1))

    while s < n
        k = k + 1
        s = s + loss_pdf(nu_l, k, xbasis)
    end
    loss_norm = loss_pdf(nu_l, k, xbasis)

    return k, loss_norm
end

#################################################

function error_prep(loss, dephase, nu_l, xbasis, sample_no)
    plus = xbasis[1]
    min = xbasis[2]

    zero = xbasis[5]
    one = xbasis[6]

    a_b = xbasis[4]
    n_b = xbasis[3]

    E(x, phi, nu) = (((1 - exp(-nu))^(x/2))/(sqrtfactorial(big(x))))*exp((-nu_l/2 + phi)*dense(n_b))*a_b^(x)

    row, col = size(loss[:, :, 1])

    err_prep_plus = fill(E(loss[1, 1, 1], dephase[1, 1, 1], nu_l)*plus, (row, col, sample_no))
    err_prep_min = fill(E(loss[1, 1, 1], dephase[1, 1, 1], nu_l)*min, (row, col, sample_no))
    err_prep_zero = fill(E(loss[1, 1, 1], dephase[1, 1, 1], nu_l)*plus, (row, col, sample_no))
    err_prep_one = fill(E(loss[1, 1, 1], dephase[1, 1, 1], nu_l)*min, (row, col, sample_no))

    for k = 1:sample_no
        for i = 1:row
            for j = 1:col
                err_prep_plus[i, j, k] = E(loss[i, j, k], dephase[i, j, k], nu_l)*plus
                err_prep_min[i, j, k] = E(loss[i, j, k], dephase[i, j, k], nu_l)*min
                err_prep_zero[i, j, k] = E(loss[i, j, k], dephase[i, j, k], nu_l)*zero
                err_prep_one[i, j, k] = E(loss[i, j, k], dephase[i, j, k], nu_l)*one
            end
        end
    end

    return err_prep_plus, err_prep_min, err_prep_zero, err_prep_one
end

####### Find exp values for easy pdf construction

function error_exp_1(err_prep_plus, err_prep_min)

    row, col, sample_no = size(err_prep_plus)
    err_exp_plus = fill(dagger(err_prep_plus[1, 1, 1])*err_prep_plus[1, 1, 1], (row, col, sample_no))
    err_exp_min = fill(dagger(err_prep_min[1, 1, 1])*err_prep_min[1, 1, 1], (row, col, sample_no))

    for k = 1:sample_no
        for i = 1:row
            for j = 1:col
                err_exp_plus[i, j, k] = dagger(err_prep_plus[i, j, k])*err_prep_plus[i, j, k]
                err_exp_min[i, j, k] = dagger(err_prep_min[i, j, k])*err_prep_min[i, j, k]
            end
        end
    end

    return err_exp_plus, err_exp_min
end

function error_exp_2(err_prep_zero, err_prep_one)

    row, col, sample_no = size(err_prep_zero)

    err_exp_zero = zeros(Complex(Float64), (row, col, sample_no))
    err_exp_one = zeros(Complex(Float64), (row, col, sample_no))
    err_exp_zo = zeros(Complex(Float64), (row, col, sample_no))
    err_exp_oz = zeros(Complex(Float64), (row, col, sample_no))

    for k = 1:sample_no
        for i = 1:row
            for j = 1:col
                err_exp_zero[i, j, k] = dagger(err_prep_zero[i, j, k])*err_prep_zero[i, j, k]
                err_exp_one[i, j, k] = dagger(err_prep_one[i, j, k])*err_prep_one[i, j, k]
                err_exp_zo[i, j, k] = dagger(err_prep_zero[i, j, k])*err_prep_one[i, j, k]
                err_exp_oz[i, j, k] = dagger(err_prep_one[i, j, k])*err_prep_zero[i, j, k]
            end
        end
    end

    return err_exp_zero, err_exp_one, err_exp_zo, err_exp_oz
end
