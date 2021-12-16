
using BesselK, SpecialFunctions, ArbNumerics
include("shared.jl")

function wbesselk(v,x)
  try
    return BesselK._besselk(v,x)
  catch
    return NaN
  end
end

rbesselk(v,x) =  Float64(ArbNumerics.besselk(ArbFloat(v), ArbFloat(x)))
abesselk(v,x) =  SpecialFunctions.besselk(v, x)

const VGRID = range(0.25, 10.0, length=101) # to avoid integer v.
const XGRID = range(0.0,  50.0, length=201)[2:end] # to avoid zero x.

const BASELINE = [rbesselk(z[1], z[2]) for z in Iterators.product(VGRID, XGRID)]
const AMOS     = [abesselk(z[1], z[2]) for z in Iterators.product(VGRID, XGRID)]
const OURSOL   = [wbesselk(z[1], z[2]) for z in Iterators.product(VGRID, XGRID)]
const TOLS_A   = atolfun.(zip(BASELINE, AMOS))
const TOLS_U   = atolfun.(zip(BASELINE, OURSOL))
const TOLS_AU  = rtolfun.(zip(AMOS,     OURSOL))

@assert iszero(length(findall(isnan, OURSOL))) "There were NaNs in our attempts!"

# quick simple test to find the worst atol when B(x) < 1:
let ix = findall(x->x<=one(x), BASELINE)
  res = findmax(abs,  TOLS_U[ix])
  _ix = res[2]
  println("Worst atol when true besselk(v,x) <= 1: $(res[1])")
  println("True value:                             $(BASELINE[ix][_ix])")
end

gnuplot_save_matrix!("../plotdata/atols_amos.csv",      TOLS_A,  VGRID, XGRID)
gnuplot_save_matrix!("../plotdata/atols_ours.csv",      TOLS_U,  VGRID, XGRID)
gnuplot_save_matrix!("../plotdata/rtols_amos_ours.csv", TOLS_AU, VGRID, XGRID)

