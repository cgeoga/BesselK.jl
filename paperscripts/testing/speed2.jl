
using LinearAlgebra, BesselK, BenchmarkTools, Printf, StaticArrays
include("shared.jl")

const GRIDN  = 24
const GRID1D = range(0.0, 1.0, length=GRIDN)
const PTS    = map(x->SVector{2,Float64}(x[1], x[2]), 
                   vec(collect.(Iterators.product(GRID1D, GRID1D))))

function checkcovmat(fn, p, v)
  _p     = pv(1.0, p, v)
  M      = assemble_matrix(fn, PTS, _p)
  em     = eigmin(M)
  Mf     = cholesky!(M, check=false)
  if issuccess(Mf)
    ld = logdet(Mf)
  else
    ld = NaN
  end
  (issuccess(Mf) ? "S" : "F", em, ld)
end

const VRANGE = (0.4, 1.25, 3.5)
const PRANGE = (0.01, 1.0, 100.0)

# For now, I might even keep the timing in seconds.
for (j, (v, p)) in enumerate(Iterators.product(VRANGE, PRANGE))
  _p  = pv(1.0, p, v)
  # test 1: assembly time.
  t_u = @belapsed assemble_matrix(matern_us,   $PTS, $_p) samples=8
  t_a = @belapsed assemble_matrix(matern_amos, $PTS, $_p) samples=8
  # test 2: Cholesky success or failure:
  (s_u, em_u, ld_u) = checkcovmat(matern_us,   p, v)
  (s_a, em_a, ld_a) = checkcovmat(matern_amos, p, v)
  # PRINTING:
  # (v,p) pair
  # times (u,p)
  # eigmins (u,p)
  # eigmin difference
  # logdets (u.p)
  # logdet difference
  if j != length(VRANGE)*length(PRANGE)
    @printf "(%1.3f, %1.2f) & %1.1e & %1.1e & %1.2e & %1.2e & %1.2e & %1.2e & %1.2e & %1.2e \\\\\n" p v t_u t_a em_u em_a abs(em_u-em_a) ld_u ld_a abs(ld_u-ld_a)
  else
    @printf "(%1.3f, %1.2f) & %1.1e & %1.1e & %1.2e & %1.2e & %1.2e & %1.2e & %1.2e & %1.2e\n" p v t_u t_a em_u em_a abs(em_u-em_a) ld_u ld_a abs(ld_u-ld_a)
  end
end

