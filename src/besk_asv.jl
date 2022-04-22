
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

