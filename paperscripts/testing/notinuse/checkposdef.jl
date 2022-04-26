
using LinearAlgebra
include("shared.jl")

const GRIDN = 24

function checkcovmat(fn, p, v; display_small_piece=false)
  grid1d = range(0.0, 1.0, length=GRIDN)
  pts    = vec(collect.(Iterators.product(grid1d, grid1d)))
  M      = assemble_matrix(fn, pts, [1.0, p, v])
  em     = eigmin(M)
  if display_small_piece # just a debugging thing, really.
    println("Principle 5x5 minor:")
    display(M[1:5, 1:5])
  end
  Mf = cholesky!(M, check=false)
  (issuccess(Mf), em)
end
translate_result(res) = res ? :SUCCESS : :FAILURE

const VRANGE = (0.4, 0.55, 0.755, 0.99, 1.01, 1.55, 2.05, 3.05, 4.05)
const PRANGE = (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)

for (v, p) in Iterators.product(VRANGE, PRANGE)
  println("\n(v,p) = ($v, $p):")
  print("AMOS:")
  amos_succ = :PLACEHOLDER
  us_succ   = :PLACEHOLDER
  (amos_em, us_em) = (0.0, 0.0)
  try
    (amos_succ, amos_em) = checkcovmat(matern_amos, p, v)
    println("$(translate_result(amos_succ)), $amos_em")
  catch
    amos_succ = :FAILURE
    println(amos_succ)
  end

  print("US:  ")
  try
    (us_succ, us_em) = checkcovmat(matern_us, p, v)
    println("$(translate_result(us_succ)), $us_em")
  catch
    us_succ = :FAILURE
    println(us_succ)
  end
  println("Eigmin difference: $(abs(amos_em-us_em))")
  println("Eigmin rtol:       $(abs(amos_em-us_em)/abs(amos_em))")

  if amos_succ && !us_succ
    println("######################")
    println("!!!!WE FAILED WHERE AMOS SUCCEEDED, INVESTIGATE!!!")
    println("######################")
  end
end

