
@inline isnearint(v, tol) = abs(v-round(v)) < tol

# TODO (cg 2021/10/22 17:04): The cutoffs chosen here are actually to some
# degree chosen as a balance between pointwise accuracy AND AD derivative
# accuracy. For example, the series gets worse faster for derivatives than it
# does for raw evals. Not entirely sure what to make of that, but here we are.
#
# TODO (cg 2021/10/27 10:28): the weak point is clearly with |z| between, say, 8
# and 15. We only get atols on the order of 1e-12 for that, which is not quite
# good enough to declare total victory. I can only manage to match the speed of
# AMOS in that part of the domain, too.
function _besselk(v, x, maxit=100, tol=1e-12, order=6)
  iszero(x) && return Inf
  (real(v) < 0) && return _besselk(-v, x, maxit, tol, order)
  # TODO (cg 2021/11/16 16:44):  this is not the right way to test if you're
  # doing AD. But I don't want to hard-bake the ForwardDiff.Dual type in here.
  is_ad = !(v isa AbstractFloat)
  # Special cases:
  #
  # Note that we do a special check for using the exponential improvement in
  # the asymptotic expansion if we're doing AD.
  half_int_branch  = isinteger(v-1/2)
  #whole_int_branch = (isinteger(v) && (abs(x) < 7.0))
  half_int_branch  && return _besselk_as(v, x, Int(ceil(v)), is_ad) # AS is exact.
  #whole_int_branch && return _besselk_temme(v, x, maxit, tol, false) 
  #
  # General cases:
  #
  # TODO (cg 2021/11/01 18:03): These branches are not perfect. If you go into
  # ./testing/accuracy.jl and track down the largest rtols between this code and
  # AMOS, you will be able to fiddle around with what version you use and get a
  # better rtol/atol. But I'm at the point where whenever I tweak something like
  # that, something else gets worse and makes it a wash. I think for the time
  # being I have to stop playing with things.
  if abs(x) < 8.5 # (x < 9)
    if (v > 2.85) || isnearint(v, 0.01) || ((x > 4) && is_ad)
      return _besselk_temme(v,x,maxit,tol,false) # direct series.
    else
      return _besselk_ser(v,x,maxit,tol,false) # direct series.
    end
  elseif abs(x) < 15.0
    return _besselk_asv(v,x,12) # uniform large order expn.
  elseif abs(x) < 30.0
    return _besselk_asv(v,x,8) 
  else
    if abs(v) > 1.5
      return _besselk_asv(v,x,8) 
    else
      return _besselk_as(v,x,order,is_ad) 
    end
  end
end

# A very simple wrapper that will use AMOS when possible, but fall back to this
# code. So now you can use AD on this but still get AMOS for direct evals.
#
# TODO (cg 2021/11/05 16:57): what's the most sensible naming thing here?
# Calling it besselk and not exporting it seems reasonable enough, but users
# will obviously want to import it. So not obvious what's best to do here.
function adbesselk(v, x, maxit=100, tol=1e-12, order=5)
  if v isa AbstractFloat
    SpecialFunctions.besselk(v, x)
  else
    _besselk(v, x, maxit, tol)
  end
end

# Not exactly a taylor series, but accurate enough.
#
# TODO (cg 2021/11/10 18:33): an enhancement here would be something that also
# worked for integers. But that seems hard. I did the whole Temme thing because
# integers are hard. But as it turns out, we really need it, because the Temme
# recursion depends on (x^v)*besselk(v,x) ->_{z->0} some finite number. But for
# v=0, which is needed to compute for _all_ integer v, that doesn't hold! An
# expansion like this that is valid for all v, including integer v, but is ALSO
# AD-compatible would take care of that entirely, because the K_{v+1}(x) term in
# the Temme series doesn't need f0.
@inline function besselkxv_t0(v, x)
  gv  = gamma(v)
  _2v = 2^(v-1)
  cof = (gv*_2v, zero(v), (_2v/4)*gv/(v-1), zero(v), (_2v/16)*gv/(v*v - 3*v + 2))
  evalpoly(x, cof)
end

function adbesselkxv(v, x, maxit=100, tol=1e-12, order=5)
  (iszero(v) && iszero(x)) && return Inf
  if !isinteger(v) && (abs(x) <= 1e-8) # use Taylor at zero.
    return besselkxv_t0(v, x)
  elseif (abs(x) < 8.5) && (v < 3.0) && !isnearint(v, 0.01)
    return _besselk_ser(v, x, maxit, tol, true)
  elseif (abs(x) < 8.5)
    return _besselk_temme(v, x, maxit, tol, true)
  else
    return adbesselk(v, x, maxit, tol, order)*(x^v)
  end
end

