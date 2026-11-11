const LARGE_SQRT_PID2 = 1.2533141373155003

# LARGE_A[i] = 1/(8i), LARGE_B[i] = (2i-1)^2/(8i), i = 1:21
const LARGE_A = (
    0.125,
    0.0625,
    0.041666666666666664,
    0.03125,
    0.025,
    0.020833333333333332,
    0.017857142857142856,
    0.015625,
    0.013888888888888888,
    0.0125,
    0.011363636363636364,
    0.010416666666666666,
    0.009615384615384616,
    0.008928571428571428,
    0.008333333333333333,
    0.0078125,
    0.007352941176470588,
    0.006944444444444444,
    0.006578947368421052,
    0.00625,
    0.005952380952380952,
)

const LARGE_B = (
    0.125,
    0.5625,
    1.0416666666666667,
    1.53125,
    2.025,
    2.5208333333333335,
    3.017857142857143,
    3.515625,
    4.013888888888889,
    4.5125,
    5.011363636363637,
    5.510416666666667,
    6.009615384615385,
    6.508928571428571,
    7.008333333333334,
    7.5078125,
    8.007352941176471,
    8.506944444444445,
    9.006578947368421,
    9.50625,
    10.005952380952381,
)

@inline function series(::Val{N}, u4::NTuple{2, T}, invx) where {N, T}
    u1, u2 = u4
    p1 = 1.0
    p2 = 1.0
    @inbounds for i in (N - 1):-1:1
        ai = LARGE_A[i]
        bi = LARGE_B[i]
        p1 = muladd(muladd(u1, ai, -bi) * invx, p1, 1.0)
        p2 = muladd(muladd(u2, ai, -bi) * invx, p2, 1.0)
    end
    return (p1, p2)
end

@inline function asexp_fixed_order(::Val{N}, a, x) where {N}
    b = a < 1.0 ? 1.0 - a : a - 1.0
    invx = 1.0 / x
    s = series(Val(N), (4.0 * a * a, 4.0 * b * b), invx)
    sc = LARGE_SQRT_PID2 * sqrt(invx)
    return (s[1] * sc, s[2] * sc)
end

@inline function besselkx_large_pair(a, x)
    return if x < 30.0
        asexp_fixed_order(Val(22), a, x)
    elseif x < 100.0
        asexp_fixed_order(Val(16), a, x)
    else
        asexp_fixed_order(Val(10), a, x)
    end
end

@inline function besselk_large_pair(a, x)
    kx_a, kx_am1 = besselkx_large_pair(a, x)
    e = exp(-x)
    return (kx_a * e, kx_am1 * e)
end
