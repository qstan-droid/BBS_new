using QuantumOptics
using Distributions

function majority(out, sample_no)
    maj_out = zeros(Bool, 1, sample_no)
    row, col = size(out)

    for i = collect(1:col)
        plus = 0
        minus = 0

        for j = collect(1:row)

            if out[j, i] == true
                plus += 1
            else
                minus += 1
            end
        end

        if plus > minus
            maj_out[i] = true
        else
            maj_out[i] = false
        end
    end

    return maj_out
end

function meas_outcome(meas, N_ord)

    phi = mod2pi(mod2pi(angle(meas)) + pi/(N_ord*2))
    if phi == 2*pi
        phi = 0
    end

    k = 0
    edge = false

    while k < N_ord

        if phi > 0 + 2*k*pi/N_ord && phi < pi/N_ord + 2*k*pi/N_ord
            return true
            edge = true
        elseif phi > pi/N_ord + 2*pi*k/N_ord && phi < (2*pi/N_ord) + 2*k*pi/N_ord
            return false
            edge = true
        end
        k += 1
    end

    if edge == false
        coin = rand(1:2)
        if coin == 1
            return true
        else
            return false
        end
    end
end

####### Optimized decoder ########
# This decoder will weight each decision on how likely it was for the beta to be wrong
# The more likely it is, the less weight the beta has

function likelihood_decoder(samples, outcome, sample_no, block_size, xbasis)

    final_out = zeros(block_size, sample_no)
    plus_cat = xbasis[1]
    min_cat = xbasis[2]
    b = xbasis[7]

    for j = collect(1:sample_no)    # cycle through
        pos_log_prob = 0.0
        neg_log_prob = 0.0

        for i = collect(1:block_size)

            if outcome[i, j] == true

                # if the outcome is positive, we check the likelihood that it was a negative
                pos_prob = (dagger(coherentstate(b, samples[i, j]))*min_cat)*(dagger(min_cat)*coherentstate(b, samples[i, j]))
                pos_log_prob = pos_log_prob - log(norm(pos_prob))
            else
                # if outcome is negative, we check the likelihood that it was positive
                neg_prob = (dagger(coherentstate(b, samples[i, j]))*plus_cat)*(dagger(plus_cat)*coherentstate(b, samples[i, j]))
                neg_log_prob = neg_log_prob - log(norm(neg_prob))
            end
        end

        # Decide on the outcome
        if pos_log_prob > neg_log_prob
            final_out[j] = true
        elseif pos_log_prob < neg_log_prob
            final_out[j] = false
        else
            coin = rand(1:2)
            if coin == 1
                final_out[j] = true
            else
                final_out[j] = false
            end
        end
    end

    return final_out
end

#### Syndrome based decoder ####

# take all samples from second rail (or first) and find the average angle that they were rotated by, then rotating back by that angle

function optimise_dephase(samples, N_ord)

    average_angle = atan(sum(sin(mod2pi(angle(samples[i]))) for i = 1:length(samples))/length(samples), sum(cos(mod2pi(angle(samples[i]))) for i = 1:length(samples))/length(samples))

    # Find the nearest angle
    gap = pi/N_ord
    diff_list = zeros(2*N_ord+1)

    for k = collect(1:2*N_ord+1)
        diff_list[k] = mod2pi(gap*(k-1)) - mod2pi(average_angle)
    end

    abs_diff_list = zeros(length(diff_list))
    for i = collect(1:length(diff_list))
        abs_diff_list[i] = abs(diff_list[i])
    end

    k_min = argmin(abs_diff_list)

    # rotate each sample back
    for i = collect(1:length(samples))

        if diff_list[k_min] < 0
            samples[i] = samples[i]*exp(-abs_diff_list[k_min]*1im)
        elseif diff_list[k_min] > 0
            samples[i] = samples[i]*exp(abs_diff_list[k_min]*1im)
        elseif diff_list[k_min] == 0
            samples[i] = samples[i]
        end

        #samples[i] = samples[i]*exp(abs_diff_list[k_min]*1im)
    end

    return samples, average_angle
end
