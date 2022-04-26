
using SpecialFunctions, ForwardDiff, FiniteDifferences
include("shared.jl")

const BIG_FD    = central_fdm(10,1)
const BIG_FD_O2 = central_fdm(10,2)

xvbesselkdv(v,x)  = ForwardDiff.derivative(_v->BesselK.adbesselkxv(_v,x), v)
xvbesselkdv2(v,x) = ForwardDiff.derivative(_v->xvbesselkdv(_v,x), v)

beskxv(v,x) = besselk(v,x)*(x^v)
dxvbesselkdv(v, x)  = BIG_FD(_v->beskxv(_v,x), v)
dxvbesselkdv2(v, x) = BIG_FD_O2(_v->beskxv(_v,x), v)

fastdxvbesselkdv2(v, x) = (beskxv(v+2e-6, x) - 2*beskxv(v+1e-6, x) + beskxv(v, x))/1e-12

const VGRID = range(0.25, 10.0, length=101)        # to avoid integer v.
const XGRID = range(0.0,  50.0, length=201)[2:end] # to avoid zero x.

const BASELINE = [dxvbesselkdv2(z[1], z[2])     for z in Iterators.product(VGRID, XGRID)]
const OURSOL   = [xvbesselkdv2(z[1], z[2])      for z in Iterators.product(VGRID, XGRID)]
const FASTFD   = [fastdxvbesselkdv2(z[1], z[2]) for z in Iterators.product(VGRID, XGRID)]
const TOLS     = atolfun.(zip(BASELINE, OURSOL))
const FDTOLS   = atolfun.(zip(BASELINE, FASTFD))
const DIFTOLS  = log10.(TOLS) .- log10.(FDTOLS)

gnuplot_save_matrix!("../plotdata/atols_deriv2_xv_fd.csv",    FDTOLS, VGRID, XGRID)
gnuplot_save_matrix!("../plotdata/atols_deriv2_xv_ad.csv",      TOLS, VGRID, XGRID)
gnuplot_save_matrix!("../plotdata/atols_deriv2_xv_fdad.csv", DIFTOLS, VGRID, XGRID)

