using QuantumOptics
using Distributions

# Prepare cat codes of order N_ord with amplitude alpha
function code_prep(N_ord, dim, alpha, code)

    if code == "cat"
        b = FockBasis(dim)

        zero_cat = normalize!(sum(coherentstate(b, alpha * exp((1im * i * pi) / N_ord)) for i = 0:2*N_ord-1))
        one_cat = normalize!(sum(coherentstate(b, alpha * exp((1im * i * pi) / N_ord)) * (-1)^i for i = 0:2*N_ord-1))

        plus_cat = (zero_cat + one_cat)/sqrt(2)
        min_cat = (zero_cat - one_cat)/sqrt(2)

        n_b = number(b)
        a_b = destroy(b)

        prep_state = [plus_cat, min_cat, n_b, a_b, zero_cat, one_cat, b, alpha, dim]

    elseif code == "binomial"
        b = FockBasis(dim)

        plus_bin = (sum(sqrt((1/(2^(alpha-1)))*(binomial(alpha, k)))*fockstate(b, k*N_ord) for k = 0:alpha))/sqrt(2)
        min_bin = (sum((-1)^k * sqrt((1/(2^(alpha-1)))*(binomial(alpha, k)))*fockstate(b, k*N_ord) for k = 0:alpha))/sqrt(2)

        zero_bin = (plus_bin + min_bin)/(sqrt(2))
        one_bin = (plus_bin - min_bin)/(sqrt(2))

        n_b = number(b)
        a_b = destroy(b)

        prep_state = [plus_bin, min_bin, n_b, a_b, zero_bin, one_bin, b, alpha, dim]
    end

    return prep_state
end

function binom_prep(N_ord, dim, N)

    b = FockBasis(dim)

    plus_bin = (sum(sqrt(binomial(N+1, i))*fockstate(b, (N_ord)*i) for i = 0:N+1)/(sqrt(2^(N+1))))
    min_bin = (sum((-1)^(i)*sqrt(binomial(N+1, i))*fockstate(b, (N_ord)*i) for i = 0:N+1)/(sqrt(2^(N+1))))

    zero_bin = (plus_bin + min_bin)/(sqrt(2))
    one_bin = (plus_bin - min_bin)/(sqrt(2))

    n_b = number(b)
    a_b = destroy(b)

    prep_state = [plus_bin, min_bin, n_b, a_b, zero_bin, one_bin, b, N]

    return prep_state
end
