using QuantumOptics
using Distributions

function loss_sample(err_place, nu_loss, m, n, xbasis, sample_no)
    loss_samples = zeros(Int64, (n, m, sample_no))
    loss_norm = fill(1.0, (n, m, sample_no))

    if err_place == true
        for k = 1:sample_no
            for i = 1:n
                for j = 1:m
                    loss_samples[i, j, k], loss_norm[i, j, k] = discrete_loss_sampling(nu_loss, xbasis)
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

function error_propagation(loss_1, dephase_1, loss_2, dephase_2, block, N_ord, sample_no)
    n = block[1]
    m = block[2]
    no_rep = block[3]

    for l = 1:sample_no
        for i = 1:n
            # do the spread from block 1 to block 2
            for j = 1:m
                for k = 1:no_rep
                    dephase_2[i, m*(k - 1)+j, l] = dephase_2[i, m*(k - 1)+j, l] - pi*loss_1[i, j, l]/(N_ord[1]*N_ord[2])
                end
            end

            # then do spread from block 2 to block 1
            for j = 1:m
                for k = 1:no_rep
                    dephase_1[i, j, l] = dephase_1[i, j, l] - pi*loss_2[i, m*(k - 1)+j, l]/(N_ord[1]*N_ord[2])
                end
            end
        end
    end

    return loss_1, dephase_1, loss_2, dephase_2
end

#### ERROR SAMPLING ####

function loss_pdf(nu_l, k_l, xbasis)
    A = (((1 - exp(-nu_l))^(k_l/2))/sqrt(factorial(big(k_l))))*exp(-nu_l*dense(xbasis[3])/2)*xbasis[4]^k_l

    #ans = (dagger(A*xbasis[5])*A*xbasis[5] + dagger(A*xbasis[6])*A*xbasis[6])/2
    ans = dagger(A*xbasis[1])*A*xbasis[1]
    
    return ans
end

function discrete_loss_sampling(nu_l, xbasis)
    # initialise
    k::Int64 = 0
    s = loss_pdf(nu_l, k, xbasis)
    n = rand(Uniform(0, 1))

    while abs(s) < n
        k = k + 1
        s = s + loss_pdf(nu_l, k, xbasis)
    end
    loss_norm = loss_pdf(nu_l, k, xbasis)

    return k, loss_norm
end

#################################################

function error_prep(loss, dephase, nu_l, xbasis, block_no)

    basis_1 = xbasis[5]
    basis_2 = xbasis[6]

    a_b = xbasis[4]
    n_b = xbasis[3]

    if block_no == 1
        E = function(x::Int64, phi, nu) 
            (((1 - exp(-nu))^(x/2))/(sqrt(factorial(big(x)))))*exp(1im*phi*dense(n_b))*exp(-nu/2*dense(n_b))*(a_b^x)
        end
    elseif block_no == 2
        E = function(x::Int64, phi, nu) 
            (((1 - exp(-nu))^(x/2))/(sqrt(factorial(big(x)))))*exp(-nu/2*dense(n_b))*(a_b^x) * exp(1im*phi*dense(n_b))
        end
    end

    row, col, sample_no = size(loss)

    err_prep_1 = fill(E(loss[1, 1, 1], dephase[1, 1, 1], nu_l)*basis_1, (row, col, sample_no))
    err_prep_2 = fill(E(loss[1, 1, 1], dephase[1, 1, 1], nu_l)*basis_2, (row, col, sample_no))

    for k = 1:sample_no
        for i = 1:row
            for j = 1:col
                err_prep_1[i, j, k] = E(loss[i, j, k], dephase[i, j, k], nu_l)*basis_1
                err_prep_2[i, j, k] = E(loss[i, j, k], dephase[i, j, k], nu_l)*basis_2
            end
        end
    end

    return err_prep_1, err_prep_2
end

####### Find exp values for easy pdf construction

function error_exp(err_prep_1, err_prep_2)

    row, col, sample_no = size(err_prep_1)

    err_exp_1 = zeros(Complex{Float64}, (row, col, sample_no))
    err_exp_2 = zeros(Complex{Float64}, (row, col, sample_no))
    err_exp_12 = zeros(Complex{Float64}, (row, col, sample_no))
    err_exp_21 = zeros(Complex{Float64}, (row, col, sample_no))

    for k = 1:sample_no
        for i = 1:row
            for j = 1:col
                err_exp_1[i, j, k] = dagger(err_prep_1[i, j, k])*err_prep_1[i, j, k]
                err_exp_2[i, j, k] = dagger(err_prep_2[i, j, k])*err_prep_2[i, j, k]
                err_exp_12[i, j, k] = dagger(err_prep_1[i, j, k])*err_prep_2[i, j, k]
                err_exp_21[i, j, k] = dagger(err_prep_2[i, j, k])*err_prep_1[i, j, k]
            end
        end
    end

    return err_exp_1, err_exp_2, err_exp_12, err_exp_21
end
