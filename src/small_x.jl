# Float64 small-argument evaluation for 0 < x <= 0.5.
#
# Reduce |v| to n + mu with |mu| <= 1/2, approximate the starting
# coefficients in u = mu^2, and recover the x-dependence in t = x^2/4.

const _SMALLX_CENTERED_P0_NUM = (
  1.1406467498033415,
  -0.017074030274400433,
  -0.0011246486847055234,
  4.9051558121754117e-5,
  1.1742234054166335e-6,
  6.583426133079264e-9,
)

const _SMALLX_CENTERED_Q0_NUM = (
  -0.7078813046532797,
  -0.0404003344955058,
  -7.39273747454712e-5,
  2.3115878439095988e-5,
  2.481298337100841e-7,
  3.757989108319689e-10,
)

const _SMALLX_CENTERED_P0Q0_DEN = (
  0.9995241497421135,
  -0.1562242539067626,
  0.000953179635782112,
  0.0001513540869946055,
  -1.9721600119915404e-6,
  9.562757160461534e-9,
)

const _SMALLX_CENTERED_SINHC_COEFFS = (
  1.0,
  0.16666666666666666,
  0.008333333333333333,
  0.0001984126984126984,
  2.7557319223985893e-6,
  2.505210838544172e-8,
  1.6059043836821613e-10,
  7.647163731819816e-13,
  2.8114572543455206e-15,
)

const _SMALLX_CENTERED_COSH_COEFFS = (
  1.0,
  0.5,
  0.041666666666666664,
  0.001388888888888889,
  2.48015873015873e-5,
  2.755731922398589e-7,
  2.08767569878681e-9,
  1.1470745597729725e-11,
  4.779477332387385e-14,
)

@inline function _smallx_centered_p0q0(mu::Float64)
  u = mu*mu
  s = muladd(8.0, u, -1.0)
  invden = inv(evalpoly(s, _SMALLX_CENTERED_P0Q0_DEN))
  p0 = evalpoly(s, _SMALLX_CENTERED_P0_NUM)*invden
  q0 = evalpoly(s, _SMALLX_CENTERED_Q0_NUM)*invden
  p0, q0, u
end

@inline function _smallx_centered_geometry(mu::Float64, x::Float64)
  L = log(2.0/x)
  a = mu*L
  ea = exp(a)

  if abs(a) < 0.5
    z = a*a
    sinh_over_mu = L*evalpoly(z, _SMALLX_CENTERED_SINHC_COEFFS)
    cosha = evalpoly(z, _SMALLX_CENTERED_COSH_COEFFS)
  else
    ema = inv(ea)
    cosha = (ea+ema)/2
    sinh_over_mu = (ea-ema)/(2mu)
  end

  sinh_over_mu, cosha, ea
end

@inline function _smallx_centered_pair_core(
  mu::Float64,
  x::Float64,
  ::Val{N},
) where N
  p, q, u = _smallx_centered_p0q0(mu)
  S, C, ea = _smallx_centered_geometry(mu, x)

  f = muladd(S, p, C*q)
  ak = ea*muladd(mu, q, p)
  Kmu = f
  xKmu1 = ak
  t = x*x/4
  tk = 1.0

  @inbounds for k in 1:N
    kf = Float64(k)
    invden = inv(kf*(kf*kf-u))
    pn = muladd(kf, p, u*q)*invden
    qn = muladd(kf, q, p)*invden

    tk *= t
    f = muladd(S, pn, C*qn)
    ak = ea*muladd(mu, qn, pn)
    Kmu = muladd(tk, f, Kmu)
    xKmu1 = muladd(tk, muladd(-2.0kf, f, ak), xKmu1)

    p = pn
    q = qn
  end

  Kmu, xKmu1/x
end

@inline function _smallx_centered_pair(mu::Float64, x::Float64)
  if x > 0.44
    return _smallx_centered_pair_core(mu, x, Val(7))
  elseif x > 0.25
    return _smallx_centered_pair_core(mu, x, Val(6))
  elseif x > 0.1
    return _smallx_centered_pair_core(mu, x, Val(5))
  elseif x > 0.03
    return _smallx_centered_pair_core(mu, x, Val(4))
  elseif x > 0.004
    return _smallx_centered_pair_core(mu, x, Val(3))
  else
    return _smallx_centered_pair_core(mu, x, Val(2))
  end
end

@inline function _besselkx_centered_geometry(v::Float64, x::Float64)
  if !(isfinite(v) && isfinite(x)) || x <= 0.0 || x > 0.5
    return Bessels.besselkx(abs(v), x)
  end

  av = abs(v)
  n = floor(Int, av+0.5)
  mu = av-n

  if mu <= -0.5
    n -= 1
    mu += 1.0
  end

  Kmu, Kmu1 = _smallx_centered_pair(mu, x)

  if n == 0
    K = Kmu
  elseif n == 1
    K = Kmu1
  else
    Km1 = Kmu
    K = Kmu1
    two_over_x = 2.0/x

    @inbounds for j in 1:(n-1)
      nu = mu+j
      Km1, K = K, muladd(nu*two_over_x, K, Km1)
    end
  end

  exp(x)*K
end

const _BESSELK_SMALL_X_CENTERED_MAX_ORDER = 24.0
const _BESSELK_SMALL_X_HALF_MAX_ORDER = 23.5
const _BESSELK_SMALL_X_CENTERED_MAX_X = 0.5

# Exact order-1/2 value followed by upward recurrence.
@inline function _besselkx_small_half_integer(
  av::Float64,
  x::Float64,
)
  n = unsafe_trunc(Int, av)
  previous = sqrt(π/(2.0*x))
  n == 0 && return previous

  inverse_x = inv(x)
  current = muladd(inverse_x, previous, previous)
  n == 1 && return current

  @inbounds for k in 1:(n-1)
    coefficient = (2.0*Float64(k)+1.0)*inverse_x
    previous, current =
      current, muladd(coefficient, current, previous)
  end

  current
end

@inline function _besselk_small_centered(
  v::Float64,
  x::Float64,
)
  av = abs(v)
  n = unsafe_trunc(Int, av)
  nf = Float64(n)

  scaled = if av == nf+0.5 &&
              av <= _BESSELK_SMALL_X_HALF_MAX_ORDER
    _besselkx_small_half_integer(av, x)
  else
    _besselkx_centered_geometry(av, x)
  end

  exp(-x)*scaled
end
