using QuantumOptics

##################################################
# Finds the coefficients behind each term
function find_coeff(meas_exp_1, meas_exp_2, block_size)

    row, col, sample_no = size(samples_1)
    rep = block_size[3]
    P = zeros(Complex{Float64}, (4, sample_no))

    for k = 1:sample_no
        A_plus = (prod(prod(sqrt(abs(meas_exp_1[1][i, j, k]*pi)) for j = 1:col) + prod(sqrt(abs(meas_exp_1[2][i, j, k]*pi)) for j = 1:col) for i = 1:row) +
                    prod(prod(sqrt(abs(meas_exp_1[1][i, j, k]*pi)) for j = 1:col) - prod(sqrt(abs(meas_exp_1[2][i, j, k]*pi)) for j = 1:col) for i = 1:row))/(sqrt(2)^(row + 1))
        A_min = (prod(prod(sqrt(abs(meas_exp_1[1][i, j, k]*pi)) for j = 1:col) + prod(sqrt(abs(meas_exp_1[2][i, j, k]*pi)) for j = 1:col) for i = 1:row) -
                    prod(prod(sqrt(abs(meas_exp_1[1][i, j, k]*pi)) for j = 1:col) - prod(sqrt(abs(meas_exp_1[2][i, j, k]*pi)) for j = 1:col) for i = 1:row))/(sqrt(2)^(row + 1))

        B_plus = (prod(prod(sqrt(abs(meas_exp_2[1][i, j, k]*pi)) for j = 1:col) + prod(sqrt(abs(meas_exp_2[2][i, j, k]*pi)) for j = 1:col) for i = 1:row*rep))/(sqrt(2)^(row*rep))
        B_min = (prod(prod(sqrt(abs(meas_exp_2[1][i, j, k]*pi)) for j = 1:col) - prod(sqrt(abs(meas_exp_2[2][i, j, k]*pi)) for j = 1:col) for i = 1:row*rep))/(sqrt(2)^(row*rep))

        P[1, k] = A_plus*B_plus
        P[2, k] = A_plus*B_min
        P[3, k] = A_min*B_plus
        P[4, k] = A_min*B_min
    end

    return P
end

##################################################
# Finds the average fidelity

function fid_ave(outcomes_1, outcomes_2, P)

    # initialise the state
    b = SpinBasis(1//2)
    I_d = identityoperator(b)
    X_d = sigmax(b)
    Z_d = sigmaz(b)

    zero = spinup(b)
    one = spindown(b)

    psi_ini = (tensor(zero, zero) + tensor(one, one))/sqrt(2)
    psi_ini_dm = dm(psi_ini)

    sample_no = size(outcomes_1)
    fid_list = zeros(Float64, sample_no)

    for k = 1:sample_no
        psi_out = normalize((P[1, k]*psi_ini + P[2, k]*tensor(I_d, X_d)*psi_ini + P[3, k]*tensor(I_d, Z_d)*psi_ini + P[4, k]*tensor(I_d, X_d*Z_d)*psi_ini)/2)

        if outcomes_1[k] == true && outcomes_2[k] == true
            correction = tensor(I_d, I_d)
        elseif outcomes_1[k] == true && outcomes_2[k] == false
            correction = tensor(I_d, X_d)
        elseif outcomes_1[k] == false && outcomes_2[k] == true
            correction = tensor(I_d, Z_d)
        elseif outcomes_1[k] == false && outcomes_2[k] == false
            correction = tensor(I_d, Z_d*X_d)
        end

        psi_corr = correction*psi_out
        psi_corr_dm = dm(psi_corr)

        # record fidelity
        fid_list[k] = real(fidelity(psi_ini_dm, psi_corr_dm))^2
    end

    # calculate average fidelity
    fid_ave = sum(fid_list)/sample_no
    println("average_fid: ", fid_ave)
    return ave_fid
end
##################################################
