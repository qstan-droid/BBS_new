using QuantumOptics
using Distributions
using BenchmarkTools
include("measurement.jl")
include("errors.jl")

function decoding(samples_1, samples_2, N_ord, block_size, err_place, err_info, sample_no, decode_type, measure, bias, xbasis, code, H_mn)
    # initialise empty arrays
    outcomes_1 = zeros(Bool, sample_no)
    outcomes_2 = zeros(Bool, sample_no)

    # how to decode?
    if decode_type == "naive"

        # decode each qubit individually
        samples_out_1 = naive_decode(samples_1, N_ord, block_size, measure[1], bias[1], decode_type, err_info, 1, xbasis)
        samples_out_2 = naive_decode(samples_2, N_ord, block_size, measure[2], bias[2], decode_type, err_info, 2, xbasis)

        # decide outcome through these
        outcomes_1 = block_decode(samples_out_1, 1, block_size)
        outcomes_2 = block_decode(samples_out_2, 2, block_size)
    elseif decode_type == "naive_ave_bias"

        # decode each qubit individually
        samples_out_1 = naive_decode(samples_1, N_ord, block_size, measure[1], bias[1], decode_type, err_info, 1, xbasis)
        samples_out_2 = naive_decode(samples_2, N_ord, block_size, measure[2], bias[2], decode_type, err_info, 2, xbasis)

        # decide outcome through these
        outcomes_1 = block_decode(samples_out_1, 1, block_size)
        outcomes_2 = block_decode(samples_out_2, 2, block_size)
    elseif decode_type == "mlnc"

        # take the two samples and some information about the system and
        outcomes_1 = max_like_no_corr(samples_1, N_ord, err_info, xbasis, measure[1], code[1], 1, bias)
        outcomes_2 = max_like_no_corr(samples_2, N_ord, err_info, xbasis, measure[2], code[2], 2, bias)
    elseif decode_type == "ml"

        # Take the two samples and find the maximum likelihood
        outcomes_1, outcomes_2 = max_like_decoder(samples_1, samples_2, N_ord, err_info, xbasis, measure, H_mn)
        #outcomes_1, outcomes_2 = ml_decoder_new(samples_1, samples_2, N_ord, err_info, xbasis, measure, code)
    elseif decode_type == "ml_ave"

        # like maximum likelihood but we average out the correlations between the measurements (should be faster?)
        #outcomes_1 = max_like_decoder_ave(samples_1, N_ord, err_info, xbasis, measure, 1)
        #outcomes_2 = max_like_decoder_ave(samples_2, N_ord, err_info, xbasis, measure, 2)

        # find the outcomes 
        meas_decode_1 = ml_ave(samples_1, N_ord, err_info, xbasis, measure, 1, block_size, H_mn)
        meas_decode_2 = ml_ave(samples_2, N_ord, err_info, xbasis, measure, 2, block_size, H_mn)

        # block_decode the outcomes
        outcomes_1 = block_decode(meas_decode_1, 1, block_size)
        outcomes_2 = block_decode(meas_decode_2, 2, block_size)
    end

    return outcomes_1, outcomes_2
end

#########################################################
# different functions, different decoding types

function naive_decode(samples, N_ord, block_size, measure, bias, decode_type, err_info, block_no, xbasis)
    # decide the outcome for each qubit individually
    row, col, sample_no = size(samples)
    samples_out = zeros(Int64, (row, col, sample_no))
    if decode_type == "naive_ave_bias"
        # bias is -phi as we want to add the bias angle to whatever result we get
        if block_no == 1
            bias = -find_ave_angle(err_info[3], N_ord, xbasis[2])
        elseif block_no == 2
            bias = -find_ave_angle(err_info[1], N_ord, xbasis[1])
        end
    end

    for k = 1:sample_no
        for i = 1:row
            for j = 1:col
                samples_out[i, j, k] = meas_outcome(samples[i, j, k], N_ord[block_no], measure, bias)
            end
        end
    end

    return samples_out
end

function max_like_decoder(samples_1, samples_2, N_ord, err_info, xbasis, measure, H_mn)
    part = zeros(4)
    row, col, sample_no = size(samples_1)

    xbasis_1 = xbasis[1]
    xbasis_2 = xbasis[2]

    N_ord_1 = N_ord[1]
    N_ord_2 = N_ord[2]

    # prepare measurement operators
    meas_op_1 = measurement_operator(measure[1], xbasis_1, N_ord_1, H_mn)
    meas_op_2 = measurement_operator(measure[2], xbasis_2, N_ord_2, H_mn)
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
            (((1-exp(-nu_loss_1))^(p1/2))/sqrt(factorial(big(p1))))*exp(-(nu_loss_1/2 + 1im*(p2*pi/(N_ord_1*N_ord_2)))*dense(n_b_1))*a_b_1^p1
            
        end
        println("A1")
    else
        A = function(p1, p2) 
            (((1-exp(-nu_loss_1))^(p1/2))/sqrt(factorial(big(p1))))*exp(-nu_loss_1*dense(n_b_1)/2)*a_b_1^p1
        end
        println("A2")
    end

    if nu_loss_1 != 0
        println("B1")
        B = function(p1, p2) 
            #(((1-exp(-nu_loss_2))^(p1/2))/sqrt(factorial(big(p1))))*exp(-nu_loss_2*dense(n_b_2)/2)*a_b_2^p1 * exp(-1im*(p2*pi/(N_ord_1*N_ord_2))*dense(n_b_2))
            exp(-1im*(p2*pi/(N_ord_1*N_ord_2))*dense(n_b_2))
        end
    else
        println("B2")
        B = function(p1, p2)
            (((1-exp(-nu_loss_2))^(p1/2))/sqrt(factorial(big(p1))))*exp(-nu_loss_2*dense(n_b_2)/2)*a_b_2^p1
        end
    end

    # prepare the superoperator version of the loss channel applied to all the states
    #sigma = [tensor(plus_1, plus_2), tensor(plus_1, min_2), tensor(min_1, plus_2), tensor(min_1, min_2)]
    #l_max = 50
    #sigma_prep = []

    #for i = 1:4
    #    push!(sigma_prep, sum(spre for l = 0:l_max))
    #end

    # find all likelihoods and choose the max

    # prepare the errored state
    lmax = findmin([xbasis_1[8]*N_ord[1], xbasis_2[8]*N_ord[2]])[1] + 2 # channel limit

    for i = collect(1:sample_no)
        
        meas_ops_1 = meas_ops[1](samples_1[i])
        meas_ops_2 = meas_ops[2](samples_2[i])
        #meas_ops_comb = tensor(meas_ops_1, meas_ops_2)

        if nu_loss_1 == 0.0 && nu_loss_2 == 0.0
            part[1] = norm(tr(meas_ops_1*plus_1)*tr(meas_ops_2*plus_2))
            part[2] = norm(tr(meas_ops_1*plus_1)*tr(meas_ops_2*min_2))
            part[3] = norm(tr(meas_ops_1*min_1)*tr(meas_ops_2*plus_2))
            part[4] = norm(tr(meas_ops_1*min_1)*tr(meas_ops_2*min_2))
        else
            A_prep_plus = [] 
            A_prep_min = []
            B_prep_plus = []
            B_prep_min = []
    
            for j = 0:lmax
                push!(A_prep_plus, tr(meas_ops_1*A(j, 0)*plus_1*dagger(A(j, 0)))) 
                push!(A_prep_min, tr(meas_ops_1*A(j, 0)*min_1*dagger(A(j, 0))))
    
                push!(B_prep_plus, tr(meas_ops_2*B(0, j)*plus_2*dagger(B(0, j)))) 
                push!(B_prep_min, tr(meas_ops_2*B(0, j)*min_2*dagger(B(0, j))))
            end

            part[1] = norm(sum(A_prep_plus[j]*B_prep_plus[j] for j = 1:lmax+1))
            part[2] = norm(sum(A_prep_plus[j]*B_prep_min[j] for j = 1:lmax+1))
            part[3] = norm(sum(A_prep_min[j]*B_prep_plus[j] for j = 1:lmax+1))
            part[4] = norm(sum(A_prep_min[j]*B_prep_min[j] for j = 1:lmax+1))
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
# new faster max_likelihood decoder

function ml_decoder_new(samples_1, samples_2, N_ord, err_info, xbasis, measure, code)
    row, col, sample_no = size(samples_1)

    # outcomes for majority
    maj_1 = zeros(Bool, (row, col, sample_no))
    maj_2 = zeros(Bool, (row, col, sample_no))

    # parts
    xbasis_1 = xbasis[1]
    xbasis_2 = xbasis[2]

    parts = zeros(Float64, 4)
    l_max = findmin([xbasis_1[8]*N_ord[1], xbasis_2[8]*N_ord[2]])[1] + 2
    phi = 0

    for ind = 1:sample_no

        if err_info[1] != 0.0 && err_info[3] == 0.0
            parts[1] = norm(sum(A_trace(samples_1[ind], measure[1], code[1], err_info[1], N_ord, xbasis_1, "plus", k, 0, phi)*B_trace(samples_2[ind], measure[2], code[2], err_info[3], N_ord, xbasis_2, "plus", 0, k, phi) for k = 0:l_max))
            parts[2] = norm(sum(A_trace(samples_1[ind], measure[1], code[1], err_info[1], N_ord, xbasis_1, "plus", k, 0, phi)*B_trace(samples_2[ind], measure[2], code[2], err_info[3], N_ord, xbasis_2, "min", 0, k, phi) for k = 0:l_max))
            parts[3] = norm(sum(A_trace(samples_1[ind], measure[1], code[1], err_info[1], N_ord, xbasis_1, "min", k, 0, phi)*B_trace(samples_2[ind], measure[2], code[2], err_info[3], N_ord, xbasis_2, "plus", 0, k, phi) for k = 0:l_max))
            parts[4] = norm(sum(A_trace(samples_1[ind], measure[1], code[1], err_info[1], N_ord, xbasis_1, "min", k, 0, phi)*B_trace(samples_2[ind], measure[2], code[2], err_info[3], N_ord, xbasis_2, "min", 0, k, phi) for k = 0:l_max))
        elseif err_info[1] == 0.0 && err_info[3] != 0.0
            parts[1] = norm(sum(A_trace(samples_1[ind], measure[1], code[1], err_info[1], N_ord, xbasis_1, "plus", 0, k, phi)*B_trace(samples_2[ind], measure[2], code[2], err_info[3], N_ord, xbasis_2, "plus", k, 0, phi) for k = 0:l_max))
            parts[2] = norm(sum(A_trace(samples_1[ind], measure[1], code[1], err_info[1], N_ord, xbasis_1, "plus", 0, k, phi)*B_trace(samples_2[ind], measure[2], code[2], err_info[3], N_ord, xbasis_2, "min", k, 0, phi) for k = 0:l_max))
            parts[3] = norm(sum(A_trace(samples_1[ind], measure[1], code[1], err_info[1], N_ord, xbasis_1, "min", 0, k, phi)*B_trace(samples_2[ind], measure[2], code[2], err_info[3], N_ord, xbasis_2, "plus", k, 0, phi) for k = 0:l_max))
            parts[4] = norm(sum(A_trace(samples_1[ind], measure[1], code[1], err_info[1], N_ord, xbasis_1, "min", 0, k, phi)*B_trace(samples_2[ind], measure[2], code[2], err_info[3], N_ord, xbasis_2, "min", k, 0, phi) for k = 0:l_max))
        elseif err_info[1] != 0.0 && err_info[3] != 0.0
            parts[1] = norm(sum(sum(A_trace(samples_1[ind], measure[1], code[1], err_info[1], N_ord, xbasis_1, "plus", k, j, phi)*B_trace(samples_2[ind], measure[2], code[2], err_info[3], N_ord, xbasis_2, "plus", j, k, phi) for k = 0:l_max) for j = 0:l_max))
            parts[2] = norm(sum(sum(A_trace(samples_1[ind], measure[1], code[1], err_info[1], N_ord, xbasis_1, "plus", k, j, phi)*B_trace(samples_2[ind], measure[2], code[2], err_info[3], N_ord, xbasis_2, "min", j, k, phi) for k = 0:l_max) for j = 0:l_max))
            parts[3] = norm(sum(sum(A_trace(samples_1[ind], measure[1], code[1], err_info[1], N_ord, xbasis_1, "min", k, j, phi)*B_trace(samples_2[ind], measure[2], code[2], err_info[3], N_ord, xbasis_2, "plus", j, k, phi) for k = 0:l_max) for j = 0:l_max))
            parts[4] = norm(sum(sum(A_trace(samples_1[ind], measure[1], code[1], err_info[1], N_ord, xbasis_1, "min", k, j, phi)*B_trace(samples_2[ind], measure[2], code[2], err_info[3], N_ord, xbasis_2, "min", j, k, phi) for k = 0:l_max) for j = 0:l_max))
        elseif err_info[1] == 0.0 && err_info[3] == 0.0
            parts[1] = norm(A_trace(samples_1[ind], measure[1], code[1], err_info[1], N_ord, xbasis_1, "plus", 0, 0, phi)*B_trace(samples_2[ind], measure[2], code[2], err_info[3], N_ord, xbasis_2, "plus", 0, 0, phi))
            parts[2] = norm(A_trace(samples_1[ind], measure[1], code[1], err_info[1], N_ord, xbasis_1, "plus", 0, 0, phi)*B_trace(samples_2[ind], measure[2], code[2], err_info[3], N_ord, xbasis_2, "min", 0, 0, phi))
            parts[3] = norm(A_trace(samples_1[ind], measure[1], code[1], err_info[1], N_ord, xbasis_1, "min", 0, 0, phi)*B_trace(samples_2[ind], measure[2], code[2], err_info[3], N_ord, xbasis_2, "plus", 0, 0, phi))
            parts[4] = norm(A_trace(samples_1[ind], measure[1], code[1], err_info[1], N_ord, xbasis_1, "min", 0, 0, phi)*B_trace(samples_2[ind], measure[2], code[2], err_info[3], N_ord, xbasis_2, "min", 0, 0, phi))
        end
        max_index = findmax(parts)[2]

        if max_index == 1
            maj_1[ind] = true
            maj_2[ind] = true
        elseif max_index == 2
            maj_1[ind] = true
            maj_2[ind] = false
        elseif max_index == 3
            maj_1[ind] = false
            maj_2[ind] = true
        elseif max_index == 4
            maj_1[ind] = false
            maj_2[ind] = false
        end
    end

    return maj_1, maj_2
end

function A_trace(x, measure, code, err_info, N_ord, xbasis, pm, l1, l2, phi)
    N_ord_1 = N_ord[1]
    N_ord_2 = N_ord[2]

    nu = err_info
    K_max = xbasis[8]

    if measure == "heterodyne"
        if code == "binomial"
            # constants
            C_l = (1 - exp(-nu))^(l1/2) /sqrt(factorial(big(l1)))
            G = function(t)
                return sqrt(binomial(K_max, t)/(2^K_max))
            end
            L = function(t)
                if t*N_ord_1 - l1 >= 0
                    sqrt(factorial(big(t*N_ord_1)))/factorial(big(t*N_ord_1 - l1))
                else
                    0
                end
            end

            exp_1 = function(s, t) 
                exp(-nu*(N_ord_1*(s+t) - 2*l1)/2)
            end
            exp_2 = function(s, t)
                exp(1im*N_ord_1*(t-s)*(pi*l2/(N_ord_1*N_ord_2) - phi))
            end

            if pm == "plus"
                # ans = (C_l^2*exp(-abs(x)^2)/pi)*sum(sum(G(j)*G(k)*L(k)*L(j)*exp_1(k, j)*exp_2(k, j)*x^(j*N_ord_1 - l1)*conj(x)^(k*N_ord_1 - l1) for k = 0:K_max) for j = 0:K_max)
                ans = (C_l^2*exp(-abs(x)^2)/pi)*sum(sum(G(j)*G(k)*L(k)*L(j)*exp_1(k, j)*exp_2(k, j)*x^(j*N_ord_1 - l1) *conj(x)^(k*N_ord_1 - l1) for k = 0:K_max) for j = 0:K_max)
            elseif pm == "min"
                ans = (C_l^2*exp(-abs(x)^2)/pi)*sum(sum((-1)^(k+j) *G(j)*G(k)*L(k)*L(j)*exp_1(k, j)*exp_2(k, j)*x^(j*N_ord_1 - l1)*conj(x)^(k*N_ord_1 - l1) for k = 0:K_max) for j = 0:K_max)
            end
        elseif code == "cat"

        end

    elseif measure == "opt_phase"
        if code == "binomial"
            # constants
            C_l = (1 - exp(-nu))^(l1/2) /sqrt(factorial(big(l1)))
            G = function(t) 
                sqrt(binomial(K_max, t)/(2^K_max))
            end
            L = function(t)
                if t*N_ord_1 - l1 >= 0
                    sqrt(factorial(big(t*N_ord_1))/factorial(big(t*N_ord_1 - l1)))
                else
                    0
                end
            end
            exp_1 = function(s, t) 
                exp(-nu*(N_ord_1*(s + t) - 2*l1)/2)
            end
            exp_2 = function(s, t) 
                exp(1im*N_ord_1*(t-s)*(pi*l2/(N_ord_1*N_ord_2) - phi))
            end

            if pm == "plus"
                ans = (C_l^2/(K_max*N_ord_1 + 1))*sum(sum(G(j)*G(k)*L(k)*L(j)*exp_1(k, j)*exp_2(k, j)*exp(1im*x*(j*N_ord_1 - l1))*exp(-1im*x*(k*N_ord_1 - l1)) for k = 0:K_max) for j = 0:K_max)
            elseif pm == "min"
                ans = (C_l^2/(K_max*N_ord_1 + 1))*sum(sum((-1)^(j + k) *G(j)*G(k)*L(k)*L(j)*exp_1(k, j)*exp_2(k, j)*exp(1im*x*(j*N_ord_1 - l1))*exp(-1im*x*(k*N_ord_1 - l1)) for k = 0:K_max) for j = 0:K_max)
            end
        elseif code == "cat"

        end
    end

    return ans
end

function B_trace(x, measure, code, err_info, N_ord, xbasis, pm, l1, l2, phi)
    # a lazy swap when we switch from block 1 to block 2
    N_ord_1 = N_ord[2]
    N_ord_2 = N_ord[1]

    nu = err_info
    K_max = xbasis[8]

    if measure == "heterodyne"
        if code == "binomial"
            # constants
            C_l = (1 - exp(-nu))^(l1/2) /sqrt(factorial(big(l1)))
            G = function(t) 
                sqrt(binomial(K_max, t)/(2^K_max))
            end
            L = function(t)
                if t*N_ord_1 - l1 >= 0
                    sqrt(factorial(big(t*N_ord_1)))/factorial(big(t*N_ord_1 - l1))
                else
                    0
                end
            end
            exp_1 = function(s, t) 
                exp(-nu*(N_ord_1*(s+t) - 2*l1)/2)
            end
            exp_2 = function(s, t)
                exp(1im*N_ord_1*(t-s)*(pi*l2/(N_ord_1*N_ord_2) - phi))
            end

            if pm == "plus"
                ans = (C_l^2*exp(-abs(x)^2)/pi)*sum(sum(G(j)*G(k)*L(k)*L(j)*exp_1(k, j)*exp_2(k, j)*x^(j*N_ord_1 - l1)*conj(x)^(k*N_ord_1 - l1) for k = 0:K_max) for j = 0:K_max)
            elseif pm == "min"
                ans = (C_l^2*exp(-abs(x)^2)/pi)*sum(sum((-1)^(k+j) *G(j)*G(k)*L(k)*L(j)*exp_1(k, j)*exp_2(k, j)*x^(j*N_ord_1 - l1)*conj(x)^(k*N_ord_1 - l1) for k = 0:K_max) for j = 0:K_max)
            end
        elseif code == "cat" 

        end

    elseif measure == "opt_phase"
        if code == "binomial"
            # constants
            C_l = (1 - exp(-nu))^(l1/2) /sqrt(factorial(big(l1)))
            G = function(t) 
                sqrt(binomial(K_max, t)/(2^K_max))
            end
            L = function(t)
                if t*N_ord_1 - l1 >= 0
                    sqrt(factorial(big(t*N_ord_1))/factorial(big(t*N_ord_1 - l1)))
                else
                    0
                end
            end
            exp_1 = function(s, t) 
                exp(-nu*(N_ord_1*(s + t) - 2*l1)/2)
            end
            exp_2 = function(s, t) 
                exp(1im*N_ord_1*(t-s)*(pi*l2/(N_ord_1*N_ord_2) - phi))
            end
            if pm == "plus"
                ans = (C_l^2/(K_max*N_ord_1 + 1))*sum(sum(G(j)*G(k)*L(k)*L(j)*exp_1(k, j)*exp_2(k, j)*exp(1im*x*(j*N_ord_1 - l1))*exp(-1im*x*(k*N_ord_1 - l1)) for k = 0:K_max) for j = 0:K_max)
            elseif pm == "min"
                ans = (C_l^2/(K_max*N_ord_1 + 1))*sum(sum((-1)^(k+j) *G(j)*G(k)*L(k)*L(j)*exp_1(k, j)*exp_2(k, j)*exp(1im*x*(j*N_ord_1 - l1))*exp(-1im*x*(k*N_ord_1 - l1)) for k = 0:K_max) for j = 0:K_max)
            end
        elseif code == "cat"

        end
    end

    return ans
end

#########################################################
# averaged_correlation decoder

function max_like_decoder_ave(samples, N_ord, err_info, xbasis, measure, block_no)
    part = zeros(2)
    row, col, sample_no = size(samples)

    # outcomes for majority
    maj = zeros(Bool, (row, col, sample_no))

    # N_ord unpacking
    N_ord_1 = N_ord[1]
    N_ord_2 = N_ord[2]

    # Fock space operators
    xbasis_1 = xbasis[1]
    xbasis_2 = xbasis[2]

    n_b_1 = xbasis_1[3]
    a_b_1 = xbasis_1[4]

    n_b_2 = xbasis_2[3]
    a_b_2 = xbasis_2[4]

    # prepare the states
    plus_1 = tensor(xbasis_1[1], dagger(xbasis_1[1]))
    min_1 = tensor(xbasis_1[2], dagger(xbasis_1[2]))

    plus_2 = tensor(xbasis_2[1], dagger(xbasis_2[1]))
    min_2 = tensor(xbasis_2[2], dagger(xbasis_2[2]))

    # unpack error information
    nu_loss_1 = err_info[1]
    nu_dephase_1 = err_info[2]
    nu_loss_2 = err_info[3]
    nu_dephase_2 = err_info[4]

    # prepare Kraus operators
    if nu_loss_2 != 0
        A = function(p1, p2)
            (((1-exp(-nu_loss_1))^(p1/2))/sqrt(factorial(big(p1))))*exp(-nu_loss_1*dense(n_b_1)/2)*a_b_1^p1 * exp(-1im*(p2*rep*pi/(N_ord_1*N_ord_2))*dense(n_b_2))
        end
    else
        A = function(p1, p2) 
            (((1-exp(-nu_loss_1))^(p1/2))/sqrt(factorial(big(p1))))*exp(-nu_loss_1*dense(n_b_1)/2)*a_b_1^p1
        end
    end

    if nu_loss_1 != 0
        B = function(p1, p2) 
            (((1-exp(-nu_loss_2))^(p1/2))/sqrt(factorial(big(p1))))*exp(-1im*(p2*pi/(N_ord_1*N_ord_2))*dense(n_b_2))*exp(-nu_loss_2*dense(n_b_2)/2)*a_b_2^p1 
            #exp(-1im*(p2*pi/(N_ord_1*N_ord_2))*dense(n_b_2))
        end
    else
        B = function(p1, p2)
            (((1-exp(-nu_loss_2))^(p1/2))/sqrt(factorial(big(p1))))*exp(-nu_loss_2*dense(n_b_2)/2)*a_b_2^p1
        end
    end

    # prepare the state
    lmax = findmin([xbasis_1[8]*N_ord[1], xbasis_2[8]*N_ord[2]])[1] + 2
    if block_no == 1
        meas_op = measurement_operator(measure[1], xbasis[1], N_ord[1])
        rho_A = (plus_2 + min_2)/sqrt(2)
        rho_dash_plus = ptrace(sum(sum(tensor(A(i, j)*plus_1*dagger(A(i, j)), B(j, i)*rho_A*dagger(B(j, i)))  for j = 0:lmax) for i = 0:lmax), 2)
        rho_dash_min = ptrace(sum(sum(tensor(A(i, j)*min_1*dagger(A(i, j)), B(j, i)*rho_A*dagger(B(j, i)))  for j = 0:lmax) for i = 0:lmax), 2)
    elseif block_no == 2
        meas_op = measurement_operator(measure[2], xbasis[2], N_ord[2])
        rho_A = (plus_1 + min_1)/sqrt(2)
        rho_dash_plus = ptrace(sum(sum(tensor(A(i, j)*rho_A*dagger(A(i, j)), B(j, i)*plus_2*dagger(B(j, i)))  for j = 0:lmax) for i = 0:lmax), 1)
        rho_dash_min = ptrace(sum(sum(tensor(A(i, j)*rho_A*dagger(A(i, j)), B(j, i)*min_2*dagger(B(j, i)))  for j = 0:lmax) for i = 0:lmax), 1)
    end

    for i = 1:sample_no
        meas_ops = meas_op(samples[i])
        part[1] = norm(tr(meas_ops*rho_dash_plus))
        part[2] = norm(tr(meas_ops*rho_dash_min))

        max_index = findmax(part)[2]

        if max_index == 1
            maj[i] = true
        elseif max_index == 2
            maj[i] = false
        end
    end

    return maj


end

function ml_ave(samples, N_ord, err_info, xbasis, measure, block_no, block_size, H_mn)
    row, col, sample_no = size(samples)
    outcomes = zeros(Int64, (row, col, sample_no))

    # prepare the vector for the likelihoods
    parts = zeros(Complex{Float64}, 2)

    # prepare the xbasis stuff
    xbasis_1 = xbasis[1]
    xbasis_2 = xbasis[2]

    # prepare the N_ords
    N_ord_1 = N_ord[1]
    N_ord_2 = N_ord[2]

    # prepare the rotation operator which error rate are we interested in
    rep = block_size[3]
    if block_no == 1  # interested in block 2 loss rate
        C = function(k)
            exp(-1im*k*rep*pi*dense(xbasis_1[3])/(N_ord_1*N_ord_2))
        end
        nu_l = err_info[3]
        xbasis_l = xbasis[2]
        
    elseif block_no == 2  # interested in block 1 loss rate
        C = function(k)
            exp(-1im*k*pi*dense(xbasis_2[3])/(N_ord_1*N_ord_2))
        end
        nu_l = err_info[1]
        xbasis_l = xbasis[1]
    end

    # find the maximum l_max
    l_max = findmin([xbasis_1[8]*N_ord[1], xbasis_2[8]*N_ord[2]])[1] + 2

    # prepare measurement operator and plus and min states
    meas_op = measurement_operator(measure[block_no], xbasis[block_no], N_ord[block_no], H_mn)
    plus = tensor(xbasis[block_no][1], dagger(xbasis[block_no][1]))
    min = tensor(xbasis[block_no][2], dagger(xbasis[block_no][2]))

    # prepare the states
    plus_err_state = [C(l)*plus*dagger(C(l)) for l = 0:l_max]
    min_err_state = [C(l)*min*dagger(C(l)) for l = 0:l_max]

    for n = 1:sample_no
        for i = 1:row
            for j = 1:col
                parts[1] = sum(tr(sparse(sparse(abs(loss_pdf(nu_l, l, xbasis_l))*meas_op(samples[i, j, n])*plus_err_state[l+1]))) for l = 0:l_max)
                parts[2] = sum(tr(sparse(abs(loss_pdf(nu_l, l, xbasis_l))*meas_op(samples[i, j, n])*min_err_state[l+1])) for l = 0:l_max)

                max_like_index = findmax(abs.(parts))[2]

                if max_like_index == 1
                    outcomes[i, j, n] = 1
                elseif max_like_index == 2
                    outcomes[i, j, n] = -1
                end
            end
        end
    end

    return outcomes
end

#########################################################
# averaged_correlation decoder

function max_like_no_corr(samples, N_ord, err_info, xbasis, measure, code, block_no, bias)
    row, col, sample_no = size(samples)

    # prepare an errored state channel
    parts = zeros(Complex{Float64}, 2)

    # xbasis sorting
    if block_no == 1
        xbasis_send = xbasis[2]
        xbasis_stay = xbasis[1]
    elseif block_no == 2
        xbasis_send = xbasis[1]
        xbasis_stay = xbasis[2]
    end

    # err_info, which mode are we talking about
    if block_no == 1
        err_info_this = err_info[1]
        err_info_away = err_info[3]
    elseif block_no == 2
        err_info_this = err_info[3]
        err_info_away = err_info[1]
    end
    
    # angle to include
    phi = find_ave_angle(err_info_away, N_ord, xbasis_send) 
    #phi = bias[block_no]

    l_max = xbasis_stay[8]*N_ord[block_no] + 2
    maj = zeros(Bool, (row, col, sample_no))
    
    for i = 1:sample_no
        parts[1] = sum(A_trace(samples[i], measure, code, err_info_this, [N_ord[block_no], N_ord[block_no]], xbasis_stay, "plus", j, 0, phi) for j = 0:l_max)
        parts[2] = sum(A_trace(samples[i], measure, code, err_info_this, [N_ord[block_no], N_ord[block_no]], xbasis_stay, "min", j, 0, phi) for j = 0:l_max)
    
        max_index = findmax(abs.(parts))[2]
        if max_index == 1
            maj[i] = true
        else
            maj[i] = false
        end
    end

    return maj
end

function find_ave_angle(nu, N_ord, xbasis)
    N_ord_1 = N_ord[1]
    N_ord_2 = N_ord[2]

    l_max = 50
    k_tilde = sum(norm(loss_pdf(nu, j, xbasis))*j for j = 0:l_max)

    return -pi*k_tilde/(N_ord_1*N_ord_2)
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
    maj = zeros(Bool, sample_no)

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
