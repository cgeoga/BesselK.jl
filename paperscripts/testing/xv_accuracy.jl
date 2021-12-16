
using BesselK, SpecialFunctions, ArbNumerics
include("shared.jl")

function wbesselkxv(v,x)
  try
    return BesselK.adbesselkxv(v,x)
  catch
    return NaN
  end
end

function rbesselkxv(v,x) 
  _v = ArbReal(v)
  _x = ArbReal(x)
  xv = _x^_v
  Float64(xv*ArbNumerics.besselk(ArbFloat(v), ArbFloat(x)))
end
abesselkxv(v,x) = (x^v)*SpecialFunctions.besselk(v, x)

const VGRID = range(0.25, 10.0, length=100) 
const XGRID = range(0.0,   8.0, length=201)[2:end] # since other impls can't do x=0.

const BASELINE = [rbesselkxv(z[1], z[2]) for z in Iterators.product(VGRID, XGRID)]
const AMOS     = [abesselkxv(z[1], z[2]) for z in Iterators.product(VGRID, XGRID)]
const OURSOL   = [wbesselkxv(z[1], z[2]) for z in Iterators.product(VGRID, XGRID)]
const TOLS_A   = atolfun.(zip(BASELINE, AMOS))
const TOLS_U   = atolfun.(zip(BASELINE, OURSOL))
const TOLS_AU  = rtolfun.(zip(AMOS,     OURSOL))

#=
gnuplot_save_matrix!("../plotdata/atols_amos.csv",      TOLS_A,  VGRID, XGRID)
gnuplot_save_matrix!("../plotdata/atols_ours.csv",      TOLS_U,  VGRID, XGRID)
gnuplot_save_matrix!("../plotdata/rtols_amos_ours.csv", TOLS_AU, VGRID, XGRID)
=#
