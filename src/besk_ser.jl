
function _besselk_ser(v, x, maxit, tol, modify)
  T    = promote_type(typeof(x), typeof(v))
  out  = zero(T)
  oneT = one(T)
  twoT = oneT+oneT
  # precompute a handful of things:
  xd2      = x/twoT
  xd22     = xd2*xd2
  half     = oneT/twoT
  if modify
    e2v      = exp2(v)
    xd2_v    = (x^(2*v))/e2v
    xd2_nv   = e2v
  else
    lxd2     = log(xd2)
    xd2_v    = exp(v*lxd2)
    xd2_nv   = exp(-v*lxd2)
  end
  gam_v    = gamma(v)
  gam_nv   = gamma(-v)
  gam_1mv  = -gam_nv*v # == gamma(one(T)-v)
  gam_1mnv = gam_v*v   # == gamma(one(T)+v)
  xd2_pow  = oneT
  fact_k   = oneT
  floatk   = convert(T, 0)
  (gpv, gmv) = (gam_1mnv, gam_1mv)
  # One final re-compression of a few things:
  _t1 = gam_v*xd2_nv*gam_1mv
  _t2 = gam_nv*xd2_v*gam_1mnv
  # now the loop using Oana's series expansion, with term function manually
  # inlined for max speed:
  for _j in 0:maxit
    t1   = half*xd2_pow
    tmp  = _t1/(gmv*fact_k)
    tmp += _t2/(gpv*fact_k)
    term = t1*tmp
    out += term
    ((abs(term) < tol) && _j>5) && return out
    # Use the trick that gamma(1+k+1+v) == gamma(1+k+v)*(1+k+v) to skip gamma calls:
    (gpv, gmv) = (gpv*(oneT+v+floatk), gmv*(oneT-v+floatk)) 
    xd2_pow *= xd22
    fact_k  *= (floatk+1)
    floatk  += 1.0
  end
  throw(error("$maxit iterations reached without achieving atol $tol."))
end

