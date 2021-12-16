
using BenchmarkTools, SpecialFunctions, FiniteDifferences, ForwardDiff, StaticArrays
include("shared.jl")

const FDFAST(fn, x) = (fn(x+h)-fn(x))/h
const FD2 = central_fdm(2, 1)

dbesselk_fdfast(v, x) = FDFAST(_v->besselk(_v, x), v)
dbesselk_fd2(v, x)    = FD2(_v->besselk(_v, x), v)
dbesselk_ad(v, x)     = ForwardDiff.derivative(_v->_besselk(_v, x), v)

#=
# And note zero allocations for the AD version. Pretty good. And faster by
# significant margins---like a factor of five---than even the most reckless FD.
for (v, x) in Iterators.product((0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.25, 3.0, 4.25),
                                 (1e-4, 0.25, 2.5, 5.0, 7.5, 12.0, 25.0))
  println("(v,x) = ($v, $x):")
  print("FDFAST:")
  @btime dbesselk_fdfast($v, $x)
  print("FD2:   ")
  @btime dbesselk_fd2($v, $x)
  print("AD:    ")
  @btime dbesselk_ad($v, $x)
  print("\n\n")
end
=#

# matern function:
const oo = @SVector ones(2)
const zz = @SVector zeros(2)

@inline pv(v) = @SVector [1.1, 1.1, v]
function matern_d3_ad(v)
  ForwardDiff.derivative(_v->BesselK.matern(oo, zz, pv(_v)), v)
end

function matern_d3_fd(v)
  FDFAST(_v->BesselK.matern(oo, zz, pv(_v)), v)
end
