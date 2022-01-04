
using SpecialFunctions, BesselK, ForwardDiff, StaticArrays
import ForwardDiff: derivative

# Convenience tool that will remove allocs:
@inline pv(scale, range, v) = @SVector [scale, range, v]

# Raw Bessel functions. Note also that there is the `adbesselkxv` method in
# BesselK.jl that gives you (x^v)*besselk(v,x) directly, sometimes with at least
# a slight gain to accuracy and speed. That isn't bound here, but obviously you
# could just slightly tweak this version.
R_besselk(v, x)    = BesselK.adbesselk(v, x)
R_besselk_dx(v, x) = derivative(_x->R_besselk(v, _x), x)
R_besselk_dv(v, x) = derivative(_v->R_besselk(_v, x), v)
R_besselk_dx_dx(v, x) = derivative(_x->R_besselk_dx(v, _x), x)
R_besselk_dx_dv(v, x) = derivative(_x->R_besselk_dv(v, _x), x)
R_besselk_dv_dv(v, x) = derivative(_v->R_besselk_dv(_v, x), v)

# Unlike the Julia-native version, this takes a pre-computed distance instead of
# two coordinates x and y, because if we have to take them as straight
# Vector{Float64} items then that will make allocations in the kernel calls and
# totally kill performance. In Julia this is avoided by using StaticArrays, but
# I don't think I can ask R users to deal with that interface. If you care
# enough about performance and need the additional flexibility that that
# provides, you could extend this code...or just switch to Julia more fully.
#
# Note here that the parameter order is:
# (sigma (scale), rho (range), nu (smoothness))
#
# You can off course change all of this up, but then be careful to also update
# the book-keeping stuff in the derivative functions.
#
# TODO (cg 2022/01/02 17:54): what is the most common signature for Julia stuff
# here? As it stands, the parameters can be passed in as a tuple, but should the
# R version take all parameters as scalars to again avoid using straight
# Vector{Float64}?
function matern(dist, params)
  (sg, rho, nu) = params
  iszero(dist) && return sg*sg
  arg = sqrt(2*nu)*dist/rho
  (sg*sg*(2^(1-nu))/gamma(nu))*BesselK.adbesselkxv(nu, arg)
end

# First derivatives:
matern_d1(dist, p) = matern(dist, (sqrt(2*p[1]), p[2], p[3]))
matern_d2(dist, p) = derivative(_p->matern(dist,pv(p[1],   _p, p[3])), p[2])
matern_d3(dist, p) = derivative(_p->matern(dist,pv(p[1], p[2],   _p)), p[3])

# Second derivatives:
matern_d1_d1(dist, p) = derivative(_p->matern_d1(dist, pv(_p,   p[2], p[3])), p[1])
matern_d1_d2(dist, p) = derivative(_p->matern_d1(dist, pv(p[1],   _p, p[3])), p[2])
matern_d1_d3(dist, p) = derivative(_p->matern_d1(dist, pv(p[1], p[2],   _p)), p[3])
matern_d2_d2(dist, p) = derivative(_p->matern_d2(dist, pv(p[1],   _p, p[3])), p[2])
matern_d2_d3(dist, p) = derivative(_p->matern_d2(dist, pv(p[1], p[2],   _p)), p[3])
matern_d3_d3(dist, p) = derivative(_p->matern_d3(dist, pv(p[1], p[2],   _p)), p[3])

