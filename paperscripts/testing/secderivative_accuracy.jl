
using SpecialFunctions, ForwardDiff, FiniteDifferences
include("shared.jl")

const BIG_FD    = central_fdm(10,1)
const BIG_FD_O2 = central_fdm(10,2)

besselkdv(v,x)  = ForwardDiff.derivative(_v->BesselK._besselk(_v,x), v)
besselkdv2(v,x) = ForwardDiff.derivative(_v->besselkdv(_v,x), v)

dbesselkdv(v, x)  = BIG_FD(_v->besselk(_v,x), v)
dbesselkdv2(v, x) = BIG_FD_O2(_v->besselk(_v,x), v)

fastdbesselkdv2(v, x) = (besselk(v+2e-6, x) - 2*besselk(v+1e-6, x) + besselk(v, x))/1e-12

const VGRID = range(0.25, 10.0, length=101)        # to avoid integer v.
const XGRID = range(0.0,  50.0, length=201)[2:end] # to avoid zero x.

const BASELINE = [dbesselkdv2(z[1], z[2])  for z in Iterators.product(VGRID, XGRID)]
const OURSOL   = [besselkdv2(z[1], z[2]) for z in Iterators.product(VGRID, XGRID)]
const FASTFD   = [fastdbesselkdv2(z[1], z[2])  for z in Iterators.product(VGRID, XGRID)]
const TOLS     = atolfun.(zip(BASELINE, OURSOL))
const FDTOLS   = atolfun.(zip(BASELINE, FASTFD))
const DIFTOLS  = log10.(TOLS) .- log10.(FDTOLS)

#=
let res = findmax(DIFTOLS)
  println("Worst case derivative difference: $(res[1])")
  println("Value of dfun: $(BASELINE[res[2]])")
end
=#

gnuplot_save_matrix!("../plotdata/atols_deriv2_fd.csv",    FDTOLS, VGRID, XGRID)
gnuplot_save_matrix!("../plotdata/atols_deriv2_ad.csv",      TOLS, VGRID, XGRID)
gnuplot_save_matrix!("../plotdata/atols_deriv2_fdad.csv", DIFTOLS, VGRID, XGRID)

