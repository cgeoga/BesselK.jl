
using SpecialFunctions, ForwardDiff, FiniteDifferences
include("shared.jl")

const BIG_FD = central_fdm(10,1)

dwbesselk(v,x) = ForwardDiff.derivative(_v->BesselK._besselk(_v,x,100,1e-12,false), v)

dbesselk(v, x) = BIG_FD(_v->besselk(_v,x), v)

const _h = 1e-6
fastfdbesselk(v, x) = (besselk(v+_h, x) - besselk(v, x))/_h

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

