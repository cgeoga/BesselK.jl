
# This is an amalgam of my original asymptotic series expansion and the
# improvements provided by Michael Helton and Oscar Smith in Bessels.jl, where
# this is more or less a pending PR (#48).

# To be replaced with Bessels.SQRT_PID2 when PR is merged.
const SQRT_PID2 = sqrt(pi/2)

# For now, no exponential improvement. It requires the exponential integral
# function, which would either need to be lifted from SpecialFunctions.jl or
# re-implemented. And with an order of, like, 10, this seems to be pretty
# accurate and still faster than the uniform asymptotic expansion.
function _besselk_as(v::V, x::T, order) where {V,T}
  fv = 4*v*v
  _z = x
  ser_v = one(T)
  floatj = one(T)
  ak_numv = fv   - floatj
  factj = one(T)
  twofloatj = one(T)
  eightj = T(8)
  for _ in 1:order
    # add to the series:
    term_v = ak_numv/(factj*_z*eightj)
    ser_v += term_v
    # update ak and _z:
    floatj += one(T)
    twofloatj += T(2)
    factj  *= floatj
    fourfloatj = twofloatj*twofloatj
    ak_numv *= (fv   - fourfloatj)
    _z *= x
    eightj *= T(8)
  end
  pre_multiply = SQRT_PID2*exp(-x)/sqrt(x)
  pre_multiply*ser_v
end

