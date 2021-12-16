
using SpecialFunctions, ForwardDiff, FiniteDifferences
include("shared.jl")

const BIG_FD = central_fdm(10,1)

function dwbesselk(v,x)
  try
    return ForwardDiff.derivative(_v->_besselk(_v,x), v)
  catch
    return NaN
  end
end

function dbesselk(v, x)
  BIG_FD(_v->besselk(_v,x), v)
end

const _h = 1e-6
function fastfdbesselk(v, x)
  (besselk(v+_h, x) - besselk(v, x))/_h
end

const VGRID = range(0.25, 10.0, length=101)        # to avoid integer v.
const XGRID = range(0.0,  50.0, length=201)[2:end] # to avoid zero x.

const BASELINE = [dbesselk(z[1], z[2])  for z in Iterators.product(VGRID, XGRID)]
const OURSOL   = [dwbesselk(z[1], z[2]) for z in Iterators.product(VGRID, XGRID)]
const FASTFD   = [fastfdbesselk(z[1], z[2])  for z in Iterators.product(VGRID, XGRID)]
const TOLS     = atolfun.(zip(BASELINE, OURSOL))
const FDTOLS   = atolfun.(zip(BASELINE, FASTFD))
const DIFTOLS  = log10.(TOLS) .- log10.(FDTOLS)

let res = findmax(DIFTOLS)
  println("Worst case derivative difference: $(res[1])")
  println("Value of dfun: $(BASELINE[res[2]])")
end

gnuplot_save_matrix!("../plotdata/atols_deriv_fd.csv",    FDTOLS, VGRID, XGRID)
gnuplot_save_matrix!("../plotdata/atols_deriv_ad.csv",      TOLS, VGRID, XGRID)
gnuplot_save_matrix!("../plotdata/atols_deriv_fdad.csv", DIFTOLS, VGRID, XGRID)

