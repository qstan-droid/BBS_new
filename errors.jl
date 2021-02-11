using Distributions

function loss_sample(err_place, nu_loss, m, n, xbasis)
    loss_samples = zeros(n, m)
    loss_norm = fill(1.0, (n, m))

    if err_place == true
        for i = 1:n
            for j = 1:m
                loss_samples[i, j], loss_norm[i, j] = loss_sample(nu_loss, xbasis)
            end
        end
    end

    return loss_samples, loss_norm
end

function dephase_sample(err_place, nu_dephase, m, n)
    dephase_samples = zeros(n, m)
    dephase_norm = fill(1.0, (n, m))

    if err_place == true
        for i = 1:n
            for j = 1:m
                dephase_angle =
                dephase_samples[i, j], dephase_norm[i, j] = dephase_sample(nu_dephase, xbasis)
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

###############################################3
