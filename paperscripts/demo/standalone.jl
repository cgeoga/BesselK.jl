
# Bring a few packages into scope:
using BesselK, ForwardDiff, StaticArrays
import BesselK: gamma, adbesselkxv

# Matern covariance function, using the rescaled (x^v)*besselk(v, x):
function matern(x, y, params)
  (sg, rho, nu) = params
  dist = norm(x-y)
  iszero(dist) && return sg*sg
  arg = sqrt(2*nu)*dist/rho
  (sg*sg*(2^(1-nu))/gamma(nu))*adbesselkxv(nu, arg)
end

# Slightly fancy matrix assembly and factorization:
function assemble_matrix(points, params)
  buf = Matrix{eltype(params)}(undef, length(points), length(points))
  Threads.@threads for k in eachindex(points) # nice multi-threading
    ptk = points[k]
    buf[k,k] = matern(ptk, ptk, params)
    @inbounds for j in 1:(k-1) # turns  off array bounds checking
      buf[j,k] = matern(points[j], ptk, params)
    end
  end
  cholesky!(Symmetric(buf)) # in-place Cholesky factorization to avoid heap allocs.
end

# Negative log-likelihood, with a slight trick for the quadratic form to just
# have to solve with the triangular factor once:
function nll(points, data, params)
  Sigma = assemble_matrix(points, params)
  (logdet(Sigma) + sum(abs2, Sigma.U'\data))/2
end

# Sample locations and data, using stack-allocated arrays via StaticArrays.jl to
# make sure that the autodiff derivatives don't make any heap allocations.
# This is just a random example to demonstrate how to create the single-arg
# closure that you can pass to ForwardDiff.
const LOCS = [@SVector rand(2) for _ in 1:1000]
const DAT  = randn(length(LOCS))

# create single-argument closure for the log-likelihood:
objective(p) = nll(LOCS, DAT, p)

# AD-generated gradient and hessian. It's that easy! Plug in to your favorite
# optimizer and you're good to go.
objective_grad(p) = ForwardDiff.gradient(objective, p)
objective_hess(p) = ForwardDiff.hessian(objective, p)


