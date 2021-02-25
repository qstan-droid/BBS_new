include("code_prep.jl")
include("errors.jl")
include("measurement.jl")
include("decode.jl")
include("fidelity.jl")

# N_ord, dim and alpha are lists
function circuit(code, N_ord, dim, alpha, block_size, err_place, err_info, measure, decode_type, sample_no, bias)

    # begin code preparation
    xbasis_1 = code_prep(N_ord[1], dim[1], alpha[1], code[1])
    xbasis_2 = code_prep(N_ord[2], dim[2], alpha[2], code[2])

    # sample errors
    # block_size = [no_rows, no_cols, no_repetitions]
    # error_placement = [loss on block_1, dephase on block_1, loss on block_2, dephase on block_2]
    # error_info = [nu_loss_1, nu_dephase_1, nu_loss_2, nu_dephase_2]

    loss_1, loss_norm_1 = loss_sample(err_place[1], err_info[1], block_size[2], block_size[1], xbasis_1, sample_no)
    dephase_1, dephase_norm_1 = dephase_sample(err_place[2], err_info[2], block_size[2], block_size[1], sample_no)

    loss_2, loss_norm_2 = loss_sample(err_place[3], err_info[3], block_size[2]*block_size[3], block_size[1], xbasis_2, sample_no)
    dephase_2, dephase_norm_2 = dephase_sample(err_place[4], err_info[4], block_size[2]*block_size[3], block_size[1], sample_no)

    # propagate errors
    # we assume the brooks-preskill code
    loss_1, dephase_1, loss_2, dephase_2 = error_propagation(loss_1, dephase_1, loss_2, dephase_2, block_size, N_ord, sample_no)

    # do measurements on both blocks
    # measure = [measure_type_1, measure_type_2]

    # First we pre-prepare the measurements to make it faster
    # first define the Krause operator for error

    err_prep_1_plus, err_prep_1_min = error_prep(loss_1, dephase_1, err_info[1], xbasis_1, 1)
    err_prep_2_zero, err_prep_2_one = error_prep(loss_2, dephase_2, err_info[3], xbasis_2, 2)

    err_prep_1 = [err_prep_1_plus, err_prep_1_min]
    err_prep_2 = [err_prep_2_zero, err_prep_2_one]

    err_exp_1_plus, err_exp_1_min, err_exp_1_pm, err_exp_1_mp = error_exp(err_prep_1_plus, err_prep_1_min)
    err_exp_2_zero, err_exp_2_one, err_exp_2_zo, err_exp_2_oz = error_exp(err_prep_2_zero, err_prep_2_one)

    err_exp_1 = [err_exp_1_plus, err_exp_1_min, err_exp_1_pm, err_exp_1_mp]
    err_exp_2 = [err_exp_2_zero, err_exp_2_one, err_exp_2_zo, err_exp_2_oz]

    ###################################################################################################
    # generate measurement samples
    samples_1 = 0
    norms_1 = 0     # this is just to fill in
    meas_exp_1 = 0

    samples_1, norms_1, meas_exp_plus_1, meas_exp_min_1, meas_exp_pm_1, meas_exp_mp_1 = measurement_samples(err_prep_1, err_prep_2, err_exp_1, err_exp_2, 1, measure, meas_exp_1, [xbasis_1, xbasis_2], N_ord, samples_1, norms_1, code, block_size)
    meas_exp_1 = [meas_exp_plus_1, meas_exp_min_1, meas_exp_pm_1, meas_exp_mp_1]
    samples_2, norms_2, meas_exp_zero_2, meas_exp_one_2, meas_exp_zo_2, meas_exp_oz_2 = measurement_samples(err_prep_1, err_prep_2, err_exp_1, err_exp_2, 2, measure, meas_exp_1, [xbasis_1, xbasis_2], N_ord, samples_1, norms_1, code, block_size)
    meas_exp_2 = [meas_exp_zero_2, meas_exp_one_2, meas_exp_zo_2, meas_exp_oz_2]

    ###################################################################################################
    # Thus comes decoding
    # outcomes are the decoded outcomes for each block, should be an array of length(sample_no) for each
    # bias = [bias_amplitude_1, bias_amplitude_2]
    outcomes_1, outcomes_2 = decoding(samples_1, samples_2, N_ord, block_size, err_place, err_info, sample_no, decode_type, bias)

    ###################################################################################################
    # First we find the
    P = find_coeff(meas_exp_1, meas_exp_2, block_size)
    ave_fidelity = fid_ave(outcomes_1, outcomes_2, P)

    return ave_fidelity, fid_list, samples_1, samples_2
end
