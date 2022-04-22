
using LinearAlgebra
include("shared.jl")

const GRIDPTS = let GRIDN = 24
  grid1d = range(0.0, 1.0, length=GRIDN)
  pts    = vec(collect.(Iterators.product(grid1d, grid1d)))
end

const RANDPTS = [rand(2) for _ in 1:(24*24)]

function assemblecovmat(fn, p, pts)
  M      = assemble_matrix(fn, pts, p, v)
  em     = eigmin(M)
  if display_small_piece # just a debugging thing, really.
    println("Principle 5x5 minor:")
    display(M[1:5, 1:5])
  end
  Mf = cholesky(M, check=false)
  (issuccess(Mf), em)
end
translate_result(res) = res ? :SUCCESS : :FAILURE

const VRANGE = (0.4, 0.55, 0.755, 0.99, 1.01, 1.55, 2.05, 3.05, 4.05)
const PRANGE = (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)

# Cases to look more into:
# ((RAND,GRID), v=(3.05, 4.05), p=0.1)
#
# otherwise, things look good: rtols can be larger than you might expect at
# times, but that appears to be happening when both eigenvalues are exact to
# more or less eps() precision where atol is more relevant anyway.
#
for (case, v, p) in Iterators.product((:GRID, :RAND), VRANGE, PRANGE)
  pts = (case == :GRID) ? GRIDPTS : RANDPTS
  M_amos   = assemble_matrix(matern_amos, pts, (1.0, p, v))
  M_us     = assemble_matrix(matern_us,   pts, (1.0, p, v))
  # pointwise checks:
  println("($case, v=$v, p=$p):")
  println("atol difference: $(maximum(abs, M_amos - M_us))")
  println("rtol difference: $(maximum(abs, tolfun.(zip(M_amos, M_us))))")
  # eigenvalue checks, assuming successful cholesky:
  cholflag = issuccess(cholesky(M_amos, check=false))
  if cholflag
    ev_amos  = eigvals(M_amos)
    ev_us    = eigvals(M_us)
    println("largest eig atol: $(maximum(abs, ev_amos-ev_us))")
    println("largest eig rtol: $(maximum(abs, tolfun.(zip(ev_amos, ev_us))))")
  else
    println("Cholesky for M_amos failed, skipping eigenvalue checks.")
  end
end

