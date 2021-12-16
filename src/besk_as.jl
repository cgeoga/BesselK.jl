
# A maximally fast upper branch incomplete gamma function when the first
# argument is a non-positive integer.
#
# (_gamma_upper_negative_integer)
@inline function _g_u_n_i(s::Int64, x, g0, expnx)
  n    = -s
  ser  = zero(x)
  _x   = one(x)
  _sgn = copy(_x)
  facn = convert(typeof(x), factorial(n))
  fac  = facn/n
  for k in 0:(n-1)
    term  = _sgn*fac*_x
    ser  += term
    _x   *= x
    _sgn *= -one(x)
    fac  /= (n-k-1)
  end
  ((expnx/_x)*ser + _sgn*g0)/facn
end

# The best speed I was able to accomplish with this thing was just to
# compartmentalize it into its own function. Definitely not ideal, but so it is.
#
# TODO (cg 2021/10/26 15:58): get rid of all the remaining factorial calls so
# that this won't literally break for order greater than, like, 19.
function exponential_improvement(v, x, l, m)
  onex = one(x)
  twox = onex+onex
  fv     = 4*v*v
  ser    = zero(x)
  _z     = x
  floatj = onex
  ak_num = fv - floatj
  factj  = onex
  twofloatj = onex
  eightj  = 8
  expx    = exp(2*x)
  expnx   = exp(-2*x)
  expintx = -expinti(-2*x)
  ser     = expx*factorial(l-1)*_g_u_n_i(1-l, 2*x, expintx, expnx)/(2*pi)
  for j in 1:(m-1)
    # add to the series:
    s   = Int(l-j)
    _g  = expx*factorial(Int(s-1))*_g_u_n_i(1-s, 2*x, expintx, expnx)/(2*pi)
    term       = ak_num/(factj*_z*eightj)*_g
    ser       += term 
    # update ak and _z:
    floatj    += onex
    twofloatj += twox
    factj     *= floatj
    ak_num    *= (fv - twofloatj^2)
    _z        *= x
    eightj    *= 8
  end
  _sgn = isodd(l) ? -one(x) : one(x)
  ser*_sgn*twox*cospi(v)
end

function _besselk_as(v, x, order, use_remainder=true, modify=false)
  onex = one(x)
  twox = onex+onex
  fv     = 4*v*v
  ser    = zero(x)
  _z     = x
  ser    = onex #zero(x) #onex/_z
  floatj = onex
  ak_num = fv - floatj
  factj  = onex
  twofloatj = onex
  eightj = 8
  for j in 1:order
    # add to the series:
    term       = ak_num/(factj*_z*eightj)
    ser       += term 
    # update ak and _z:
    floatj    += onex
    twofloatj += twox
    factj     *= floatj
    ak_num    *= (fv - twofloatj^2)
    _z        *= x
    eightj    *= 8
  end
  if use_remainder
    _rem = exponential_improvement(v, x, Int(order+1), Int(order))
    ser += _rem
  end
  # if you're modifying as (x^v)*besselk(v,x), since the series part is pretty
  # stable numerically, what we want to deal with is the (x^v)*exp(-x). That's
  # the problem of potentially huge*tiny.
  if modify
    mulval = exp(v*log(x)-x)
  else
    mulval = exp(-x)
  end
  return sqrt(pi/(x*twox))*mulval*ser
end

