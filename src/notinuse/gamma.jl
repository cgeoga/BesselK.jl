
# A Julia-native gamma function. Pretty darn accurate and fast, and using this
# instead of SpecialFunctions.gamma would mean that the dependency can be
# removed. That doesn't seem worth it to me at this point in time, but putting
# this here in case.

const _G  = 607/128 
const _C1 = 0.99999999999999709182
const _C  = (57.156235665862923517,
            -59.597960355475491248,
             14.136097974741747174,
            -0.49191381609762019978,
            .33994649984811888699e-4,
            .46523628927048575665e-4,
           -.98374475304879564677e-4,
            .15808870322491248884e-3,
           -.21026444172410488319e-3,
            .21743961811521264320e-3,
           -.16431810653676389022e-3,
            .84418223983852743293e-4,
           -.26190838401581408670e-4,
            .36899182659531622704e-5)
const _Cl = length(_C)
const SQRT_2_PI = sqrt(2*pi)

# Lanczos approximation. 
@inline function _gamma_pos(_z::T) where{T}
  T_one    = one(T)
  T_half   = T_one/(T_one+T_one)
  z_half   = _z-T_half
  z_g_half = z_half+_G
  out = 0.99999999999999709182 
  zpj = _z 
  @inbounds for c in _C
    out += c/zpj 
    zpj += T_one 
  end
  out*SQRT_2_PI*exp(z_half*log(z_g_half) - z_g_half)
end

function _gamma(x)
  iszero(x) && return Inf
  (isinteger(x) && x < zero(x)) && return NaN
  x > zero(x) && return _gamma_pos(x)
  flx = floor(x)
  xf  = x - flx
  g   = _gamma_pos(xf) 
  _x  = xf-one(x)
  for _ in 1:Int(-flx)
    g  /= _x
    _x -= one(x)
  end
  g
end

