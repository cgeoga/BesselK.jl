
using BenchmarkTools, StaticArrays
include("shared.jl")

const oo = @SVector ones(2)
const zz = @SVector zeros(2)
const TESTING_PARAMS = (@SVector [1.25, 0.25, 0.5],
                        @SVector [1.25, 5.25, 0.5],
                        @SVector [1.25, 0.25, 0.75],
                        @SVector [1.25, 5.25, 0.75],
                        @SVector [1.25, 0.25, 1.0],
                        @SVector [1.25, 5.25, 1.0],
                        @SVector [1.25, 0.25, 1.75],
                        @SVector [1.25, 5.25, 1.75],
                        @SVector [1.25, 0.25, 3.0],
                        @SVector [1.25, 5.25, 3.0],
                        @SVector [1.25, 0.25, 4.75],
                        @SVector [1.25, 5.25, 4.75])

for pp in TESTING_PARAMS

  println("\n###")
  println("Parameters $pp")
  println("###\n")

  println("Fast finite diff, h=$h:")
  @btime matern_fdfast_d1($oo, $zz, $pp)
  @btime matern_fdfast_d2($oo, $zz, $pp)
  @btime matern_fdfast_d3($oo, $zz, $pp)

  println("Adaptive finite diff, order 2:")
  @btime matern_fd2_d1($oo, $zz, $pp)
  @btime matern_fd2_d2($oo, $zz, $pp)
  @btime matern_fd2_d3($oo, $zz, $pp)

  println("Adaptive finite diff, order 5:")
  @btime matern_fd5_d1($oo, $zz, $pp)
  @btime matern_fd5_d2($oo, $zz, $pp)
  @btime matern_fd5_d3($oo, $zz, $pp)

  println("Complex step, h=$ch:")
  try
    @btime matern_cstep_d1($oo, $zz, $pp)
    @btime matern_cstep_d2($oo, $zz, $pp)
    @btime matern_cstep_d3($oo, $zz, $pp)
  catch
    println("Failure/error. Probably at integer values.")
  end

  println("Autodiff:")
  @btime matern_ad_d1($oo, $zz, $pp)
  @btime matern_ad_d2($oo, $zz, $pp)
  @btime matern_ad_d3($oo, $zz, $pp)

end
