
const _GA        = MathConstants.γ       # because I don't like unicode...
const QUADGAMMA1 =  -2.40411380631918857 # ψ^(2)(1)
const HEXGAMMA1  = -24.88626612344087823 # ψ^(4)(1)

# special methods for cosh(mu(v,x))*(x^vf) and sinh(...).
@inline coshmuxv(v, x) = (exp2(v) + exp2(-v)*x^(2*v))/2
@inline sinhmuxv(v, x) = (exp2(v) - exp2(-v)*x^(2*v))/2

const TPCOEF = (1, 0, (pi^2)/6, 0, (7*pi^4)/360, 0, (31*pi^6)/15120)
@inline tp_taylor0(v) = evalpoly(v, TPCOEF) # trig part

const G1COEF = (_GA, 
                0, 
                (2*_GA^3 - _GA*pi^2 - 2*QUADGAMMA1)/12,
                0,
                (12*_GA^5 - 20*(_GA^3)*pi^2 + _GA*pi^4 - 120*(_GA^2)*QUADGAMMA1 + 20*(pi^2)*QUADGAMMA1 - 12*HEXGAMMA1)/1440)
@inline g1_taylor0(v) = -evalpoly(v, G1COEF)

const G2COEF = (1.0,
                0,
                (_GA^2 - (pi^2)/6)/2,
                0,
                (60*_GA^4 - 60*(_GA*pi)^2 + pi^4 - 240*_GA*QUADGAMMA1)/1440)
@inline g2_taylor0(v) = evalpoly(v, G2COEF)

const SHCOEF = (1, 0, 1/6, 0, 1/120, 0, 1/5040, 0, 1/362880)
@inline sh_taylor0(v) = evalpoly(v, SHCOEF) # sinh part

# TODO: there is still a problem when v is not zero but z is zero. The NaN might
# be mathematically correct, and it isn't a branch that adbesselkxv hits.
@inline function f0_expansion0(v::V, z, modify) where{V}
  _t  = tp_taylor0(v)
  _g1 = g1_taylor0(v)
  _g2 = g2_taylor0(v)
  if !modify
    mu  = v*log(2/z)
    _cm = cosh(mu)
    _sm = sinh(mu)
    _sh = sh_taylor0(mu)*log(2/z)
  else
    _cm = coshmuxv(v, z)
    if _iszero(v)
      # Because of the branching here, if v is zero, then I know that z is NOT zero. 
      mu  = v*log(2/z)
      _sh = (z^v)*sh_taylor0(mu)*log(2/z)
    else
      _sh = sinhmuxv(v, z)/v 
    end
  end
  _t*(_g1*_cm + _g2*_sh)
end

# EVEN Cheby coefs for computing the gamma pair thing when v is not near zero.
# (but |v| is <= 1/2).
const A2N = (1.843740587300906,
            -0.076852840844786,
             0.001271927136655,
            -0.000004971736704,
            -0.000000033126120,
             0.000000000242310,
            -0.000000000000170,
            -0.000000000000001
            )

# EVEN Cheby coefs for computing the gamma pair thing when v is not near zero.
# (but |v| is <= 1/2).
const A2Np1 = (-0.283876542276024,
                0.001706305071096,
                0.000076309597586,
               -0.000000865920800,
                0.000000001745136,
                0.000000000009161,
               -0.000000000000034)

# This gives g1 and g2 directly when v is not near zero and completely avoids
# calls to gamma functions.
@inline function temmegammas(v)
  twov = one(v)+one(v)
  _v   = twov*v
  (tv_even, tv_odd)   = (one(v), _v)
  (ser_even, ser_odd) = (A2N[1]/2, A2Np1[1]*tv_odd)
  @inbounds for n in 2:7
    # get next even value. Note that tv_even is now T_2(_v).
    tv_even   = twov*_v*tv_odd - tv_even
    # add to the even term series:
    ser_even += A2N[n]*tv_even
    # get the next odd term:
    tv_odd    = twov*_v*tv_even - tv_odd
    # add to the odd term series:
    ser_odd  += A2Np1[n]*tv_odd
  end
  # one more term for the evens:
  tv_even   = twov*_v*tv_odd - tv_even
  ser_even += A2N[8]*tv_even
  # now return:
  (ser_odd/v, ser_even)
end

function temme_pair(v, z, maxit, tol, modify=false)
  @assert abs(v) <= 1/2 "This internal routine is only for |v|<=1/2."
  # Some very low-level objects:
  onez = one(z)
  twoz = onez+onez
  zd2  = z/twoz
  _2dz = twoz/z
  # Creating the necessary things to get initial f, p, q:
  if abs(v) < 0.001
    g1 = g1_taylor0(v)
    g2 = g2_taylor0(v)
  else
    (g1, g2) = temmegammas(v)
  end
  (gp, gm) = (inv(-(g1*twoz*v - g2*twoz)/twoz), inv((g1*twoz*v + g2*twoz)/twoz))
  # p0 and q0 terms, branches for if we're modifying by (x^vf):
  if !modify 
    p0 = (exp(-v*log(zd2))/twoz)*gp
    q0 = (exp(v*log(zd2))/twoz)*gm
  else
    p0 = exp2(v-one(v))*gp
    q0 = exp2(-v-one(v))*(z^(2*v))*gm
  end
  # cosh and sinh terms for f0, branches for if we're modifying by (x^vf):
  if !modify
     mu = v*log(_2dz)
    _cm = cosh(mu)
    _sm = sinh(mu)
  else
    _cm = coshmuxv(v, z)
    _sm = sinhmuxv(v, z)
  end
  # One more branch for f0, which is to check if v is near zero or z is near two.
  if _iszero(z) && modify
    f0 = one(z)
  else
    if abs(v) < 0.001
      f0 = f0_expansion0(v, z, modify)
    elseif abs(z-2)<0.001
      _s = sh_taylor0(v*log(_2dz)) # manually plug in mu.
      #f0 = (v*pi/sinpi(v))*(g1*cosh(mu) + g2*log(_2dz)*_s)
      f0 = (v*pi/sinpi(v))*(g1*_cm + g2*log(_2dz)*_s)
    else
      # Temme's form is:
      #f0 = (v*pi/sinpi(v))*(g1*cosh(mu) + g2*log(_2dz)*sinh(mu)/mu)
      # But if I modify to this, I get rid of a near singularity as z->0:
      f0 = (v*pi/sinpi(v))*(g1*_cm + g2*_sm/v)
    end
  end
  (_f, _p, _q) = (f0, p0, q0)
  # a few other odds and ends for efficient looping:
  (ser_kv, ser_kvp1) = (f0, _p)
  (factk, _floatk)   = (onez, onez)
  (v2, _zd4, _z)     = (v*v, z*z/(twoz + twoz), z*z/(twoz + twoz))
  for k in 1:maxit
    _f   = (k*_f + _p + _q)/(_floatk^2 - v2)
    _p  /= (_floatk-v)
    _q  /= (_floatk+v)
    ck   =  _z/factk
    # update term for besselk(v, z).
    term_v    = ck*_f
    ser_kv   += term_v
    # update term for besselk(v+1, z).
    term_vp1  = ck*(_p - k*_f)
    ser_kvp1 += term_vp1
    if max(abs(term_v), abs(term_vp1)) < tol
      if !modify
        return (ser_kv, ser_kvp1*_2dz)
      else
        return (ser_kv, ser_kvp1*2) # note that I'm multiplying by a z, so cancel manually.
      end
    end
    _floatk += onez
    factk   *= _floatk
    _z      *= _zd4
  end
  throw(error("Term tolerance $tol not reached in $maxit iters for (v,z) = ($v, $z).")) 
end

# NOTE: In the modified scaling of (x^v)*besselk(v,x), there is a problem for
# integer v: (x^0)*besselk(0,x) is just besselk(0,x). And that is still Inf for
# x=0. Which really breaks the whole strategy of integer derivatives for the
# rescaled Bessel here. BUT: there is a workaround! You don't actually need
# besselk(0, x) for the modified recurrence, so we just throw away that value.
function _besselk_temme(v, z, maxit, tol, mod)
  @assert v > -1/2 "This routine does not presently handle the case of v > -1/2."
  _p   = floor(v)
  (v - _p > 1/2) && (_p += one(_p))
  vf   = v - _p 
  twov = one(v)+one(v)
  _v   = vf
  kvp2 = zero(v)
  (kv, kvp1) = temme_pair(_v, z, maxit, tol, mod)
  # check if any recurrence is necessary:
  v <= one(v)/twov          && return kv
  v <= (twov + one(v))/twov && return kvp1
  # if it is necessary, perform it:
  if mod
    # not necessarily the "right" way to handle this, but seems to stop the
    # propagation of NaNs in AD.
    #
    # a slightly different recurrence for (x^v)*besselk(v,x).
    if _iszero(z)
      # special case for z=0:
      for _ in 1:(Int(_p)-1)
        kvp2 = twov*(_v+one(v))*kvp1
        _v  += one(v)
        kv   = kvp1
        kvp1 = kvp2
      end
    else
      z2 = z*z
      for _ in 1:(Int(_p)-1)
        kvp2 = twov*(_v+one(v))*kvp1 + z2*kv
        _v  += one(v)
        kv   = kvp1
        kvp1 = kvp2
      end
    end
  else
    for _ in 1:(Int(_p)-1)
      kvp2 = ((twov*(_v+one(v)))/z)*kvp1 + kv
      _v  += one(v)
      kv   = kvp1
      kvp1 = kvp2
    end
  end
  kvp2
end

