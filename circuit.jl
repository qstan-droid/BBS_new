include("code_prep.jl")
include("errors.jl")

# N_ord, dim and alpha are lists
function circuit(code, N_ord, dim, alpha, block_size, err_place, err_info, measure, sample_no)

    # begin code preparation
    xbasis_1 = code_prep(N_ord[1], dim[1], alpha[1], code[1])
    xbasis_2 = code_prep(N_ord[2], dim[2], alpha[2], code[2])

    # sample errors
    # block_size = [no_rows, no_cols, no_repetitions]
    # error_placement = [loss on block_1, dephase on block_1, loss on block_2, dephase on block_2]
    # error_info = [nu_loss_1, nu_dephase_1, nu_loss_2, nu_dephase_2]

    loss_1, loss_norm_1 = loss_sample(err_place[1], err_info[1], block_size[2], block_size[1], xbasis_1, sample_no)
    dephase_1, dephase_norm_1 = dephase_sample(err_place[2], err_info[2], block_size[2], block_size[1], xbasis_2, sample_no)

    loss_2, loss_norm_2 = loss_sample(err_place[3], err_info[3], block_size[2]*block_size[3], block_size[1], sample_no)
    dephase_2, dephase_norm_2 = dephase_sample(err_place[4], err_info[4], block_size[2]*block_size[3], block_size[1], sample_no)

    # propagate errors
    # we assume the brooks-preskill code
    loss_1, dephase_1, loss_2, dephase_2 = error_propagation(loss_1, dephase_1, loss_2, dephase_2, block_size, N_ord, sample_no)

    # do measurements on both blocks
    # measure = [measure_op_1, measure_op_2]

    # First we pre-prepare the measurements to make it faster
    # first define the Krause operator for error


    err_prep_1_plus, err_prep_1_min = error_prep(loss_1, dephase_1, error_info[1], xbasis_1, sample_no)
    err_prep_2_plus, err_prep_2_min = error_prep(loss_2, dephase_2, error_info[2], xbasis_2, sample_no)
    err_prep_1 = [err_prep_1_plus, err_prep_1_min]
    err_prep_2 = [err_prep_2_plus, err_prep_2_min]


    samples_1 = measurement_samples(err_prep_1, err_prep_2, 1)
    samples_2 = measurement_samples(err_prep_1, err_prep_2, 2)

end
