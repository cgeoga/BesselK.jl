
# This code is taken from Bessels.jl (https://github.com/heltonmc/Bessels.jl),
# and is provided under the MIT license. Copyright Michael Helton and
# contributors.


#    Modified Bessel functions of the second kind of order zero and one
#                       besselk0, besselk1
#
#    Scaled modified Bessel functions of the second kind of order zero and one
#                       besselk0x, besselk1x
#
#    (Scaled) Modified Bessel functions of the second kind of order nu
#                       besselk, besselkx
#
#
#    Calculation of besselk0 is done in two branches using polynomial approximations [1]
#
#    Branch 1: x < 1.0 
#              besselk0(x) + log(x)*besseli0(x) = P7(x^2)
#                            besseli0(x) = [x/2]^2 * P6([x/2]^2) + 1
#    Branch 2: x >= 1.0
#              sqrt(x) * exp(x) * besselk0(x) = P22(1/x) / Q2(1/x)
#    where P7, P6, P22, and Q2 are 7, 6, 22, and 2 degree polynomials respectively.
#
#
#
#    Calculation of besselk1 is done in two branches using polynomial approximations [2]
#
#    Branch 1: x < 1.0 
#              besselk1(x) - log(x)*besseli1(x) - 1/x = x*P8(x^2)
#                            besseli1(x) = [x/2]^2 * (1 + 0.5 * (*x/2)^2 + (x/2)^4 * P5([x/2]^2))
#    Branch 2: x >= 1.0
#              sqrt(x) * exp(x) * besselk1(x) = P22(1/x) / Q2(1/x)
#    where P8, P5, P22, and Q2 are 8, 5, 22, and 2 degree polynomials respectively.
#
#
#    The polynomial coefficients are taken from boost math library [3].
#    Evaluation of the coefficients using Remez.jl is prohibitive due to the increase
#    precision required in ArbNumerics.jl. 
#
#    Horner's scheme is then used to evaluate all polynomials.
#
#    Calculation of besselk and besselkx can be done with recursion starting with
#    besselk0 and besselk1 and using upward recursion for any value of nu (order).
#
#                    K_{nu+1} = (2*nu/x)*K_{nu} + K_{nu-1}
#
#    When nu is large, a large amount of recurrence is necesary.
#    We consider uniform asymptotic expansion for large orders to more efficiently
#    compute besselk(nu, x) when nu is larger than 100 (Let's double check this cutoff)
#    The boundary should be carefully determined for accuracy and machine roundoff.
#    We use 10.41.4 from the Digital Library of Math Functions [5].
#    This is also 9.7.8 in Abramowitz and Stegun [6].
#    K_{nu}(nu*z) = sqrt(pi / 2nu) *exp(-nu*n)/(1+z^2)^1/4 * sum((-1^k)U_k(p) /nu^k)) for k=0 -> infty
#    The U polynomials are the most tricky. They are listed up to order 4 in Table 9.39
#    of [6]. For Float32, >=4 U polynomials are usually necessary. For Float64 values,
#    >= 8 orders are needed. However, this largely depends on the cutoff of order you need.
#    For moderatelly sized orders (nu=50), might need 11-12 orders to reach good enough accuracy
#    in double precision. 
#
#    However, calculation of these higher order U polynomials are tedious. These have been hand
#    calculated and somewhat crosschecked with symbolic math. There could be errors. They are listed
#    in src/U_polynomials.jl as a reference as higher orders are impossible to find while being needed for any meaningfully accurate calculation.
#    For large orders these formulas will converge much faster than using upward recurrence.
#
#    
# [1] "Rational Approximations for the Modified Bessel Function of the Second Kind 
#     - K0(x) for Computations with Double Precision" by Pavel Holoborodko     
# [2] "Rational Approximations for the Modified Bessel Function of the Second Kind 
#     - K1(x) for Computations with Double Precision" by Pavel Holoborodko
# [3] https://github.com/boostorg/math/tree/develop/include/boost/math/special_functions/detail
# [4] "Computation of Bessel Functions of Complex Argument and Large Order" by Donald E. Amos
#      Sandia National Laboratories
# [5] https://dlmf.nist.gov/10.41
# [6] Abramowitz, Milton, and Irene A. Stegun, eds. Handbook of mathematical functions with formulas, graphs, and mathematical tables. 
#     Vol. 55. US Government printing office, 1964.
#

const P1_k0(::Type{Float32}) = (
    -1.372508979104259711f-1, 2.622545986273687617f-1, 5.047103728247919836f-3
)
const Q1_k0(::Type{Float32}) = (
    1.000000000000000000f0, -8.928694018000029415f-2, 2.985980684180969241f-3
)
const P2_k0(::Type{Float32}) = (
    1.159315158f-1, 2.789828686f-1, 2.524902861f-2,
    8.457241514f-4, 1.530051997f-5
)
const P3_k0(::Type{Float32}) = (
    2.533141220f-1, 5.221502603f-1,
    6.380180669f-2, -5.934976547f-2
)
const Q3_k0(::Type{Float32}) = (
    1.000000000f0, 2.679722431f0,
    1.561635813f0, 1.573660661f-1
)
const P1_k0(::Type{Float64}) = (
    -1.372509002685546267e-1, 2.574916117833312855e-1,
    1.395474602146869316e-2, 5.445476986653926759e-4,
    7.125159422136622118e-6
)
const Q1_k0(::Type{Float64}) = (
    1.000000000000000000e+00, -5.458333438017788530e-02,
    1.291052816975251298e-03, -1.367653946978586591e-05
)
const P2_k0(::Type{Float64}) = (
    1.159315156584124484e-01, 2.789828789146031732e-01,
    2.524892993216121934e-02, 8.460350907213637784e-04,
    1.491471924309617534e-05, 1.627106892422088488e-07,
    1.208266102392756055e-09, 6.611686391749704310e-12
)
const P3_k0(::Type{Float64}) = (
    2.533141373155002416e-1, 3.628342133984595192e0,
    1.868441889406606057e1, 4.306243981063412784e1,
    4.424116209627428189e1, 1.562095339356220468e1,
    -1.810138978229410898e0, -1.414237994269995877e0,
    -9.369168119754924625e-2
)
const Q3_k0(::Type{Float64}) = (
    1.000000000000000000e0, 1.494194694879908328e1,
    8.265296455388554217e1, 2.162779506621866970e2,
    2.845145155184222157e2, 1.851714491916334995e2,
    5.486540717439723515e1, 6.118075837628957015e0,
    1.586261269326235053e-1
)
const Y_k0 = 1.137250900268554688


const Y_k1(::Type{Float32}) = 8.695471287f-2
const Y_k1(::Type{Float64}) = 8.69547128677368164e-2

const Y2_k1(::Type{Float32}) = 1.450342178f0
const Y2_k1(::Type{Float64}) = 1.45034217834472656

const P1_k1(::Type{Float32}) = (
    -3.621379531f-3, 7.131781976f-03, -1.535278300f-5
)
const P1_k1(::Type{Float64}) = (
    -3.62137953440350228e-3, 7.11842087490330300e-3,
    1.00302560256614306e-5, 1.77231085381040811e-6
)
const Q1_k1(::Type{Float32}) = (
    1.000000000f0, -5.173102701f-2, 9.203530671f-4
)
const Q1_k1(::Type{Float64}) = (
    1.00000000000000000e0, -4.80414794429043831e-2,
    9.85972641934416525e-4, -8.91196859397070326e-6
)
const P2_k1(::Type{Float32}) = (
    -3.079657469f-1, -8.537108913f-2,
    -4.640275408f-3, -1.156442414f-4
)
const P2_k1(::Type{Float64}) = (
    -3.07965757829206184e-1, -7.80929703673074907e-02,
    -2.70619343754051620e-3, -2.49549522229072008e-5
)
const Q2_k1(::Type{Float32}) = (
    one(Float32), zero(Float32)
)
const Q2_k1(::Type{Float64}) = (
    1.00000000000000000e0, -2.36316836412163098e-2,
    2.64524577525962719e-4, -1.49749618004162787e-6
)

const P3_k1(::Type{Float32}) = (
    -1.970280088f-1, 2.188747807f-2,
    7.270394756f-1, 2.490678196f-1
)
const P3_k1(::Type{Float64}) = (
    -1.97028041029226295e-1, -2.32408961548087617e0,
    -7.98269784507699938e0, -2.39968410774221632e0,
    3.28314043780858713e1, 5.67713761158496058e1,
    3.30907788466509823e1, 6.62582288933739787e0,
    3.08851840645286691e-1
)
const Q3_k1(::Type{Float32}) = (
    1.000000000f0, 2.274292882f0,
    9.904984851f-1, 4.585534549f-2
)
const Q3_k1(::Type{Float64}) = (
    1.00000000000000000e0, 1.41811409298826118e1,
    7.35979466317556420e1, 1.77821793937080859e2,
    2.11014501598705982e2, 1.19425262951064454e2,
    2.88448064302447607e1, 2.27912927104139732e0,
    2.50358186953478678e-2
)


"""
    besselk0(x::T) where T <: Union{Float32, Float64}

Modified Bessel function of the second kind of order zero, ``K_0(x)``.
"""
function besselk0(x::T) where T <: Union{Float32, Float64}
    x <= zero(T) && return throw(DomainError(x, "`x` must be nonnegative."))
    if x <= one(T)
        a = x * x / 4
        s = muladd(evalpoly(a, P1_k0(T)), inv(evalpoly(a, Q1_k0(T))), T(Y_k0))
        a = muladd(s, a, 1)
        return muladd(-a, log(x), evalpoly(x * x, P2_k0(T)))
    else
        s = exp(-x / 2)
        a = muladd(evalpoly(inv(x), P3_k0(T)), inv(evalpoly(inv(x), Q3_k0(T))), one(T)) * s / sqrt(x)
        return a * s
    end
end

"""
    besselk1(x::T) where T <: Union{Float32, Float64}

Modified Bessel function of the second kind of order one, ``K_1(x)``.
"""
function besselk1(x::T) where T <: Union{Float32, Float64}
  x <= zero(T) && return throw(DomainError(x, "`x` must be nonnegative."))
  if x <= one(T)
    z = x * x
    a = z / 4
    pq = muladd(evalpoly(a, P1_k1(T)), inv(evalpoly(a, Q1_k1(T))), Y_k1(T))
    pq = muladd(pq * a, a, (a / 2 + 1))
    a = pq * x / 2
    pq = muladd(evalpoly(z, P2_k1(T)) / evalpoly(z, Q2_k1(T)), x, inv(x))
    return muladd(a, log(x), pq)
  else
    s = exp(-x / 2)
    a = muladd(evalpoly(inv(x), P3_k1(T)), inv(evalpoly(inv(x), Q3_k1(T))), Y2_k1(T)) * s / sqrt(x)
    return a * s
  end
end

#=
If besselk0(x) or besselk1(0) is equal to zero
this will underflow and always return zero even if besselk(x, nu)
is larger than the smallest representing floating point value.
In other words, for large values of x and small to moderate values of nu,
this routine will underflow to zero.
For small to medium values of x and large values of nu this will overflow and return Inf.
=#
function _besselk_int(nu, x::T) where T <: Union{Float32, Float64, BigFloat}
  T == Float32 ? branch = 20 : branch = 50
  if nu < branch
    return up_recurrence(x, besselk0(x), besselk1(x), nu)[1]
  else
    return _besselk_asv(nu, x, 8)
  end
end


@inline function up_recurrence(x, k0, k1, nu)
  nu == 0 && return k0
  nu == 1 && return k1
  # this prevents us from looping through large values of 
  # nu when the loop will always return zero
  (iszero(k0) || iszero(k1)) && return zero(x) 
  k2 = k0
  x2 = 2 / x
  for n in 1:nu-1
    a = x2 * n
    k2 = muladd(a, k1, k0)
    k0 = k1
    k1 = k2
  end
  return k2, k0
end

