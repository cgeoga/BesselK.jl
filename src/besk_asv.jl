
@inline _eta(z) = sqrt(1+z*z)+log(z/(1+sqrt(1+z*z)))
@inline _p(z)   = inv(sqrt(1+z*z))

function Uk_polynomials(max_order)
  P0   = Polynomial([1.0])
  out  = [P0]
  mul_int  = Polynomial([1.0, 0.0, -5.0])
  mul_frnt = Polynomial([0.0, 0.0, 1.0, 0.0, -1.0])/2
  for j in 1:max_order
    Pjm1 = out[end]
    Pjm1_int = integrate(mul_int*Pjm1)/8
    Pjm1_drv = derivative(Pjm1)
    newP     = mul_frnt*Pjm1_drv + Pjm1_int - Pjm1_int(0.0)/8
    push!(out, newP)
  end
  out
end

# Should I compute these and hard-code them into this file?
if !(@isdefined UK_COEF_TUPLES)
  const UK_POLYS = Uk_polynomials(20)
end

function _besselk_vz_asv(v, z, maxorder, modify=false)
  ez  = _eta(z)
  pz  = _p(z)
  if modify
    # same trick: (x^v)*exp(-v*...) in the log.
    mulval = exp(v*log(z*v)-v*ez)
  else
    mulval = exp(-v*ez)
  end
  out = sqrt(pi/(2*v))*mulval/sqrt(sqrt(1+z*z))
  # series part:
  ser = zero(z)
  sgn = one(z)
  _v  = one(v)
  @inbounds for j in 1:maxorder
    ser += sgn*UK_POLYS[j](pz)/_v
    sgn *= -one(z)
    _v  *= v
  end
  out*ser
end

_besselk_asv(v, z, maxorder=10, modify=false) = _besselk_vz_asv(v, z/v, maxorder, modify)

#=
function _besselk_vz_asv_pairv(v1, v2, z1, z2, maxorder)
  (ez1, ez2)  = (_eta(z1), _eta(z2))
  (pz1, pz2)  = (_p(z1),   _p(z2))
  out1 = sqrt(pi/(2*v1))*exp(-v1*ez1)/sqrt(sqrt(1+z1*z1))
  out2 = sqrt(pi/(2*v2))*exp(-v2*ez2)/sqrt(sqrt(1+z2*z2))
  # series part:
  sgn = one(z1)
  (ser1, ser2) = (zero(z1), zero(z2))
  (_v1, _v2)   = (one(v1), one(v2))
  for (j, coeft) in pairs(UK_COEF_TUPLES)
    ser1 += sgn*evalpoly(pz1, coeft)/_v1
    ser2 += sgn*evalpoly(pz2, coeft)/_v2
    sgn  *= -one(z1)
    _v1  *= v1
    _v2  *= v2
    (j == maxorder) && break
  end
  (out1*ser1, out2*ser2)
end
=#


#=
_besselk_asv_pairv(v1, v2, z, maxorder=10) = _besselk_vz_asv_pairv(v1, v2, z/v1, z/v2, maxorder)
_besselk_asv_pair_naive(v1, v2, z) = (_besselk_asv(v1, z), _besselk_asv(v2, z))


# This function should not be called for |z|<5, say, because for large orders
# besselk(v_large, z) is huge and you'll lose digits.
function _besselk_asv_backrec(v, z, reclevel=3, maxorder=10)
  _vp2 = v+reclevel+1
  _vp1 = v+reclevel
  _v   = v+reclevel-one(v)
  (mp2, mp1) = _besselk_asv_pairv(_vp2, _vp1, z, maxorder)
  m = mp2 - (2*_vp1/z)*mp1
  for it in 1:(reclevel-1)
    _vp1 -= one(v)
    mp2   = mp1
    mp1   = m
    m     = mp2 - (2*_vp1/z)*mp1
  end
  m
end
=#
