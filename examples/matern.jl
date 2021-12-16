
using BesselK, SpecialFunctions, ForwardDiff, StaticArrays, FiniteDifferences

# for quick testing/debugging:
const _oo = @SVector ones(2)
const _zz = @SVector zeros(2)
const _pp = @SVector [1.1, 2.1, 1.3]

# Convenience tool that will remove allocs:
@inline pv(scale, range, v) = @SVector [scale, range, v]

# FD things to compare with:
const FDFAST(fn, x) = (fn(x+1e-6) - fn(x))/1e-6
const FDACC  = central_fdm(10,1)
const FD2    = central_fdm(2,1)


# Note that this function uses the specialized method for (x^v)*besselk(v,x).
# Which can actually be MORE accurate than using AMOS or whatever for very small
# x. Pretty interesting. Even for pretty large order (say, v \approx 5), it is
# pretty accurate. More importantly, it speeds up the AD a bit and makes THAT
# more accurate.
function matern(x, y, params)
  (sg, rho, nu) = params
  dist = norm(x-y)
  iszero(dist) && return sg*sg
  arg = sqrt(2*nu)*dist/rho
  (sg*sg*(2^(1-nu))/gamma(nu))*BesselK.adbesselkxv(nu, arg)
end

@inline beskxv(v, x) = BesselK.besselk(v, x)*(x^v)
function matern_amos(x, y, params)
  (sg, rho, nu) = params
  dist = norm(x-y)
  iszero(dist) && return sg*sg
  arg = sqrt(2*nu)*dist/rho
  (sg*sg*(2^(1-nu))/gamma(nu))*beskxv(nu, arg)
end

# Easily computed analytically, so skipping this.
matern_d1(x, y, params)    = matern(x, y, (sqrt(2*params[1]), params[2], params[3]))
matern_d1_d1(x, y, params) = matern(x, y, (sqrt(2), params[2], params[3]))

# Can be computed analytically, so we'll go ahead and assume a dedicated FD user
# has grinded out the analytical derivatives.
function matern_d2(x, y, params) 
  ForwardDiff.derivative(p2->matern(x, y, pv(params[1], p2, params[3])), params[2])
end
function matern_d1_d2(x, y, params) 
  ForwardDiff.derivative(p2->matern_d1(x, y, pv(params[1], p2, params[3])), params[2])
end
function matern_d2_d2(x, y, params) 
  ForwardDiff.derivative(p2->matern_d2(x, y, pv(params[1], p2, params[3])), params[2])
end

# First derivatives:
for (stem, fn) in ((:fdfast, FDFAST),
                   (:fd2,    FD2),
                   (:ad,     ForwardDiff.derivative))
  name3 = Symbol(:matern_, stem, :_d3)
  if stem == :ad
    @eval $name3(x,y,p) = $fn(_p->matern(x,y,pv(p[1],p[2],_p)), p[3])
  else
    @eval $name3(x,y,p) = $fn(_p->matern_amos(x,y,pv(p[1],p[2],_p)), p[3])
  end
end

# Second derivatives:
for (stem, fn) in ((:fdfast, FDFAST),
                   (:ad,     ForwardDiff.derivative))
  # first names:
  name1 = Symbol(:matern_d1)
  name2 = Symbol(:matern_d2)
  name3 = Symbol(:matern_, stem, :_d3)
  # second names:
  name13 = Symbol(:matern_, stem, :_d1, :_d3)
  name22 = Symbol(:matern_, stem, :_d2, :_d2)
  name23 = Symbol(:matern_, stem, :_d2, :_d3)
  name33 = Symbol(:matern_, stem, :_d3, :_d3)
  # d1d3:
  @eval $name13(x,y,p) = $fn(_p->$name1(x,y,pv(p[1],p[2],_p)), p[3])
  # d2d3:
  @eval $name23(x,y,p) = $fn(_p->$name2(x,y,pv(p[1],p[2],_p)), p[3])
  # d3d3:
  @eval $name33(x,y,p) = $fn(_p->$name3(x,y,pv(p[1],p[2],_p)), p[3])
end

const AD_DFNS  = (matern_d1, matern_d2, matern_ad_d3)
const AD_D2FNS = ((matern_d1_d1, matern_d1_d2, matern_ad_d1_d3),
                  (matern_d2_d2, matern_ad_d2_d3),
                  (matern_ad_d3_d3,))
      

const FD_DFNS  = (matern_d1, matern_d2, matern_fdfast_d3)
const FD_D2FNS = ((matern_d1_d1, matern_d1_d2, matern_fdfast_d1_d3),
                  (matern_d2_d2, matern_fdfast_d2_d3),
                  (matern_fdfast_d3_d3,))

