using QuantumOptics
using Distributions

function decoding(samples_1, samples_2, N_ord, block_size, err_place, err_info, sample_no, decode_type, measure)
    # initialise empty arrays
    outcomes_1 = zeros(Bool, sample_no)
    outcomes_2 = zeros(Bool, sample_no)

    # how to decode?
    if decode_type == "naive"
        # decode each qubit individually
        samples_out_1 = naive_decode(samples_1, N_ord[1], block_size, measure[1])
        samples_out_2 = naive_decode(samples_2, N_ord[2], block_size, measure[2])

        # decide outcome through these
        outcomes_1 = block_decode(samples_out_1, 1)
        outcomes_2 = block_decode(samples_out_2, 2)

    elseif decode_type == "bias"

    elseif decode_type == "max_likelihood"

    end

    return outcomes_1, outcomes_2
end

#########################################################
# different functions, different decoding types

function naive_decode(samples, N_ord, block_size, measure)
    # decide the outcome for each qubit individually
    row, col, sample_no = size(samples)
    samples_out = zeros(Int64, (row, col, sample_no))

    for k = 1:sample_no
        for i = 1:row
            for j = 1:col
                samples_out[i, j, k] = meas_outcome(samples[i, j, k], N_ord, measure, 0)
            end
        end
    end



    return samples_out
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

    if block_no == 1
        # we compute total parity of each column
        # take majority vote over all parities
        for k = 1:sample_no
            col_par = fill(1, col)

            for i = 1:col
                col_par[i] = prod(samples_out[j, i, k] for j = 1:no_row)
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
            row_par = fill(1, (rep, div(col, rep)))

            # take row parities
            for i = 1:row
                for j = 1:rep
                    row_par[j, i] = prod(samples_out[i, l] for l = (j-1)*div(col, rep)+1:j*div(col, rep))
                end
            end
            # majority vote over repetitions


            # majority vote over all column parities

        end
    end

    return maj
end
