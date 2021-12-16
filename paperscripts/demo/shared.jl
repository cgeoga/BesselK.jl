
# Setting: imagine that we have a bunch of iid replicates of very correlated
# smooth data. Let's fit it and look at the resulting CIs.

using LinearAlgebra, StableRNGs
include("../../examples/matern.jl")
include("gradient.jl")
include("efish.jl")
include("hessian.jl")

# A likelihood function, slightly specialized for speed:
function nll_replicates(parms, pts, datav)
  K    = Symmetric([matern(x, y, parms) for x in pts, y in pts])
  Kf   = cholesky!(K) # in-place
  out  = logdet(Kf)*length(datav)/2 # only compute logdet once.
  out += sum(z->sum(abs2, Kf.L\z), datav)/2 # all the solve terms.
  out
end

# For plotting a likelihood surface.
function profile_nll_replicates(parms, pts, datav)
  n    = length(first(datav))
  _p   = @SVector [one(eltype(parms)), parms[1], parms[2]]
  K    = Symmetric([matern(x, y, _p) for x in pts, y in pts])
  Kf   = cholesky!(K) # in-place
  out  = logdet(Kf)*length(datav)/2 # only compute logdet once.
  out += n*sum(z->log(sum(abs2, Kf.L\z)), datav)/2 # all the solve terms.
  out
end

# simulate some data, using the "let" syntax to avoid keeping global K after
# we're done:
const TRU_P = @SVector [1.5, 2.5, 1.3]
const SEED  = StableRNG(12345)
const N_REP = 10
const N_DAT = 512
const PTS   = [SVector{2,Float64}(rand(SEED, 2)...) for _ in 1:N_DAT]
const SIMS  = let K = Symmetric([matern(x, y, TRU_P) for x in PTS, y in PTS])
  Kf = cholesky(K)
  [Kf.L*randn(SEED, size(Kf, 2)) for _ in 1:N_REP]
end

_nll(p)  = nll_replicates(p, PTS, SIMS)
_nllh(p) = ForwardDiff.hessian(_nll, p)

function caching_nll(p, cache=nothing)
  if !isnothing(cache)
    push!(cache, deepcopy(p))
  end
  _nll(p)
end

const HIGHFD = central_fdm(10,1)
function high_fd_hessian(p)
  gfun = _p -> FiniteDifferences.grad(HIGHFD, _nll, _p)
  FiniteDifferences.jacobian(HIGHFD, gfun, p)[1]
end

grad_fd!(p, store) = gradient_replicates!(store, FD_DFNS, p, PTS, SIMS)
grad_ad!(p, store) = ForwardDiff.gradient!(store, _nll, p)
hessfd(p) = hessian_replicates(FD_DFNS, FD_D2FNS, p, PTS, SIMS)
fishfd(p) = efish_replicates(FD_DFNS, p, PTS, length(SIMS))
fishad(p) = efish_replicates(AD_DFNS, p, PTS, length(SIMS))
