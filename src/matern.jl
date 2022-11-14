
# A standard Matern covariance function. This is not the best parameterization,
# but it is the most popular, so here we are.

"""
matern(x, y, params)

computes σ² * \\mathcal{M}_{ν}(||x-y||/ρ), where params = (σ, ρ, ν) and \\mathcal{M} is the Matern covariance function, parameterized as

\\mathcal{M}_{ν}(t) = σ^2 2^{1 - ν} Γ(ν)^{-1} (\\sqrt{2 ν} t / ρ)^{ν} \\mathcal{K}_{ν}(\\sqrt{2 ν} t / ρ).

For more information, see Stein (1999), Interpolation of Spatial Data: Some Theory for Kriging.
"""
function matern(x, y, params)
  (sg, rho, nu) = (params[1], params[2], params[3])
  dist = norm(x-y)
  _iszero(dist) && return sg^2
  arg = sqrt(2*nu)*dist/rho
  (sg*sg*(2^(1-nu))/_gamma(nu))*adbesselkxv(nu, arg)
end

