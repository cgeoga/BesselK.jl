
include("../../examples/matern.jl") 
include("../plotting/gnuplot_utils.jl") 

const VGRID = range(0.25, 10.0, length=101) # to hit "near-integer" v, which is the hardest.
const XGRID = range(0.0,  30.0, length=201)[2:end]

# First and second derivatives with FD:
fdbesselkdv(v, x) = (besselk(v+1e-6, x) - besselk(v+1e-6, x))/1e-6
fdbesselkdvdv(v, x) = (besselk(v+2e-6, x) - 2*besselk(v+1e-6, x) + besselk(v, x))/1e-12

# First and second derivatives with AD:
adbesselkdv(v, x) = ForwardDiff.derivative(_v->BesselK._besselk(_v, x), v)
adbesselkdvdv(v, x) = ForwardDiff.derivative(_v->adbesselkdv(_v, x), v) 

function assemble_matrix(fn, pts, p)
  out = Array{Float64}(undef, length(pts), length(pts))
  Threads.@threads for k in 1:length(pts)
    out[k,k] = p[1]
    @inbounds for j in 1:(k-1)
      out[j,k] = fn(pts[j], pts[k], p)
    end
  end
  Symmetric(out)
end

atolfun(tru, est) = isnan(est) ? NaN : (isinf(tru) ? 0.0 : abs(tru-est))
atolfun(tru_est)  = atolfun(tru_est[1], tru_est[2])
rtolfun(tru, est) = atolfun(tru, est)/abs(tru)
rtolfun(tru_est)  = rtolfun(tru_est[1], tru_est[2])


