
# Credit: https://www.advanpix.com/2015/11/11/rational-approximations-for-the-modified-bessel-function-of-the-first-kind-i0-computations-double-precision/

const P16_COEF = (
  1.0000000000000000000000801e+00,
  2.4999999999999999999629693e-01,
  2.7777777777777777805664954e-02,
  1.7361111111111110294015271e-03,
  6.9444444444444568581891535e-05,
  1.9290123456788994104574754e-06,
  3.9367598891475388547279760e-08,
  6.1511873265092916275099070e-10,
  7.5940584360755226536109511e-12,
  7.5940582595094190098755663e-14,
  6.2760839879536225394314453e-16,
  4.3583591008893599099577755e-18,
  2.5791926805873898803749321e-20,
  1.3141332422663039834197910e-22,
  5.9203280572170548134753422e-25,
  2.0732014503197852176921968e-27,
  1.1497640034400735733456400e-29)

const P22_COEF = (
  3.9894228040143265335649948e-01,
  4.9867785050353992900698488e-02,
  2.8050628884163787533196746e-02,
  2.9219501690198775910219311e-02,
  4.4718622769244715693031735e-02,
  9.4085204199017869159183831e-02,
  -1.0699095472110916094973951e-01,
  2.2725199603010833194037016e+01,
  -1.0026890180180668595066918e+03,
  3.1275740782277570164423916e+04,
  -5.9355022509673600842060002e+05,
  2.6092888649549172879282592e+06,
  2.3518420447411254516178388e+08,
  -8.9270060370015930749184222e+09,
  1.8592340458074104721496236e+11,
  -2.6632742974569782078420204e+12,
  2.7752144774934763122129261e+13,
  -2.1323049786724612220362154e+14,
  1.1989242681178569338129044e+15,
  -4.8049082153027457378879746e+15,
  1.3012646806421079076251950e+16,
  -2.1363029690365351606041265e+16,
  1.6069467093441596329340754e+16)

# Can also use this to get a really fast _1F1(1/2,1,z), although didn't end up
# using that.
function besseli0(x)
  onex = one(x)
  xd2square = (x/(onex+onex))^2
  if abs(x) < 7.75
    return xd2square*evalpoly(xd2square, P16_COEF) + onex
  else
    return evalpoly(onex/x, P22_COEF)/(sqrt(x)*exp(-x))
  end
end

# Accurate enough up to |x|=10. 
#
# No longer used at present, but functional enough to be worth keeping.
function _besselk0_ser(x, maxit=50, tol=1e-12)
  one_x  = one(x)
  two_x  = one_x+one_x
  four_x = two_x+two_x
  xd4_2  = (x*x/four_x)
  out    = -(log(x/two_x) + MathConstants.γ)*besseli0(x)
  floatj = one_x
  _xd4_2 = xd4_2
  coef   = one_x
  factj  = one_x
  for j in 1:maxit
    term = coef*_xd4_2/(factj^2)
    out += term
    abs(term)<tol && return out
    floatj += one_x
    coef   += 1/floatj
    _xd4_2 *= xd4_2
    factj  *= floatj
  end
  throw(error("No convergence to term tolerance $tol in $maxit iterations."))
end
#_besselk1_ser(x) = -ForwardDiff.derivative(_besselk0_ser, x)

function _besselk_int(v, x)
  (bk0, bk1) = temme_pair(zero(x), x, 100, 1e-12)
  iszero(v) && return bk0#_besselk0_ser(x)
  isone(v)  && return bk1#_besselk1_ser(x)
  (vm1, bvm1, bvm2, new_val) = (1, bk1, bk0, (2/x)*bk1 + bk0)
  isequal(2, v) && return new_val
  for _ in 3:v
    vm1 += 1
    bvm2 = bvm1
    bvm1 = new_val
    new_val = (2*(vm1)/x)*bvm1 + bvm2
  end
  new_val
end

# A first derivative for integer orders. Not as fast as it could be, but still
# pretty fast. So for the moment, good enough.
#
# No longer used at present, but functional enough to be worth keeping.
function deriv_besselk_int(v, x) 
  iszero(v) && return zero(x)
  isone(v)  && return _besselk0_ser(x)/x 
  (bk0, bk1) = temme_pair(zero(x), x, 100, 1e-12)
  (vm1, bvm1, bvm2, new_val) = (1, bk1, bk0, (2/x)*bk1 + bk0)
  zd2   = x/2
  _zd2  = zd2^2
  flj   = 2.0
  factj = 2.0
  ser   = bk0/v + zd2*bk1/(v-1)
  for j in 2:(v-1)
    # add on to the series appropriately:
    term   = new_val*_zd2/(factj*(v-j))
    ser   += term
    # Find the next bessel_k(x), as new_val:
    vm1 += 1
    bvm2 = bvm1
    bvm1 = new_val
    new_val = (2*(vm1)/x)*bvm1 + bvm2
    # now increment a few things:
    flj   += one(x)
    factj *= flj
    _zd2  *= zd2
  end
  factj*ser/(2*_zd2) 
end

