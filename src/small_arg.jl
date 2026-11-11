# Julia port of small_arg.rs (Temme's method for K_mu, K_{mu+1} at small x).
#
# Rust const generics -> `::Val` dispatch; Rust stack arrays -> tuples.

const TEMME_G1 = (
    # Gamma1(mu) = (1/Gamma(1-mu) - 1/Gamma(1+mu)) / (2 mu)
    -0.5772156649015329,
    0.042002635034095515,
    0.042197734555520845,
    -0.007218943245911219,
    0.0002152416622963602,
    2.0134955675158432e-05,
    -1.1335045193083941e-06,
    -4.939613004327888e-09,
)

const TEMME_G2 = (
    # Gamma2(mu) = (1/Gamma(1-mu) + 1/Gamma(1+mu)) / 2
    1.0,
    -0.6558780715202527,
    0.16653861138218923,
    -0.009621971524606148,
    -0.0011651676432173535,
    0.0001280507201273341,
    -1.2525586982068558e-06,
    -2.0057050414754113e-07,
)

const TEMME_SINC = (
    # sin(pi mu) / (pi mu)
    1.0,
    -1.6449340668482264,
    0.8117424252833394,
    -0.19075182412149658,
    0.026147847805557536,
    -0.0023460808955200734,
    0.0001484278416841416,
    -6.972100313670779e-06,
    2.4503008901533465e-07,
)

const TEMME_SINHC = (
    # sinh(s)/s in s^2
    1.0,
    1.0 / 6.0,
    1.0 / 120.0,
    1.0 / 5040.0,
    1.0 / 362880.0,
    1.0 / 39916800.0,
    1.0 / 6227020800.0,
)

@inline function temme_pair(::Val{M}, mu, x) where {M}
    z = 0.25 * x * x
    d = -log(0.5 * x) # ln(2/x)
    sig = mu * d
    t = mu * mu

    g1 = evalpoly(t, TEMME_G1)
    g2 = evalpoly(t, TEMME_G2)
    sinc = evalpoly(t, TEMME_SINC)

    ep = exp(sig)    # (x/2)^(-mu)
    em = 1.0 / ep    # (x/2)^(+mu)
    ch = 0.5 * (ep + em)
    shc = abs(sig) < 0.25 ? evalpoly(sig * sig, TEMME_SINHC) : 0.5 * (ep - em) / sig

    # independent of the recurrence, so the divisions can overlap
    r = ntuple(Val(M)) do i
        kf = Float64(i)
        1.0 / (kf * kf - t)
    end

    f = (g1 * ch + g2 * shc * d) / sinc
    p = 0.5 * ep / (g2 - mu * g1) # (x/2)^(-mu) / (2 Gamma(1+mu))
    q = 0.5 * em / (g2 + mu * g1) # (x/2)^(+mu) / (2 Gamma(1-mu))

    c, s0, s1 = 1.0, f, p
    @inbounds for k in 1:M
        kf = Float64(k)
        rk = r[k]
        c *= z / kf
        f = muladd(kf, f, p + q) * rk
        p *= (kf + mu) * rk
        q *= (kf - mu) * rk
        s0 = muladd(c, f, s0)
        s1 = muladd(c, muladd(kf, -f, p), s1)
    end
    return (s0, s1 * 2.0 / x)
end

@inline temme_pair_binned(mu, x) =
    x <= 0.65 ? temme_pair(Val(8), mu, x) : temme_pair(Val(12), mu, x)

# Returns (K_a, K_{a-1}) for the e^x-unscaled functions, like the Rust version.
@inline function besselk_small_pair(a, x)
    a_ge_half = a >= 0.5
    mu = a_ge_half ? a - 1.0 : -a # K_{-s} = K_s, so |mu| = min(a, b) <= 1/2
    k_mu, k_mu1 = temme_pair_binned(mu, x)
    return a_ge_half ? (k_mu1, k_mu) : (k_mu, k_mu1)
end
