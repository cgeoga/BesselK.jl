
# Observations:
#
# By and large, we really smoke AMOS in speed. Especially for large arg.
# Assembling a reasonably large matrix (1024 x 1024) gets done in about half the
# time. Which may seem like nothing, but a factor of two never hurts...

using BenchmarkTools, SpecialFunctions, StaticArrays
include("../besk.jl")
include("shared.jl")

const NUs = (0.1, 0.25, 0.5, 0.75, 0.999, 1.0, 1.001, 1.25, 1.5, 1.75, 2.0,
             2.05, 2.75, 3.1, 3.8, 4.3, 4.9)

const TINY_X   = 1e-4
const SMALL_X  = 0.25
const MID_X    = 7.5
const LARGE_X  = 20.0
const Xs = (1e-4, 0.25, 5.0, 7.5, 8.5, 14.0, 17.0, 29.0, 31.0, 50.0)

const PTS  = [rand(3).*10.0 for _ in 1:1024]
const PRMS = @SVector [1.0, 1.0, 1.25]

print("\n\n")
println("##################")
println("HEAT ONE: pointwise timings.")
println("##################")
print("\n\n")
for (v, x) in Iterators.product(NUs, Xs)
  println("(v,x) = ($v, $x):")
  print("Timing for AMOS:")
  @btime SpecialFunctions.besselk($v, $x)
  print("Timing for us:")
  try
    @btime _besselk($v, $x)
  catch
    println("FAILURE/ERROR OUT.")
  end
  print("\n\n")
end

print("\n\n")
println("##################")
println("HEAT TWO: kernel matrix assembly.")
println("##################")
print("\n\n")
println("Timings for AMOS:")
@btime assemble_matrix(matern_amos, $PTS, $PRMS)
println("Timings for US:")
@btime assemble_matrix(matern_us, $PTS, $PRMS)

