
# As is discussed in the paper, you can see that the FD derivatives here are so
# inaccurate that Ipopt can't even converge. The FD Hessians are worse than
# garbage and completely unusable without very high-order and adaptive methods,
# which are so brutal to performance that they aren't even worth considering.

using Ipopt, Serialization

include("shared.jl")
include("gradient.jl")
include("ipopt_helpers.jl")

function fitter(objects_initval)
  ((case, gradfun, hesfun, maxiter), ini) = objects_initval
  println("Case: $case")
  cache = Vector{Float64}[]
  prob  = createProblem(3, [0.0, 0.0, 0.25], fill(1e22, (3,)), 0,
                        Float64[], Float64[], 0, div(3*4, 2),
                        _p->caching_nll(_p, cache),
                        (args...)->nothing,
                        gradfun,
                        (args...)->nothing,
                        (x,m,r,c,o,l,v)->ipopt_hessian(x,m,r,c,o,l,v,
                                                       hesfun,Function[],0))
  addOption(prob, "tol", 1e-5)
  addOption(prob, "max_iter", maxiter)
  prob.x = deepcopy(ini) # for safety to avoid weird persistent pointer games.
  @time status = solveProblem(prob)
  return (prob, status, case, deepcopy(prob.x), _nll(prob.x), hesfun(prob.x), cache, ini)
end

const cases = ((:FD_FISH, grad_fd!, fishfd, 100),
               (:AD_FISH, grad_ad!, fishad, 100),
               (:FD_HESS, grad_fd!, hessfd, 100), 
               (:AD_HESS, grad_ad!, _nllh,  100))
const inits = (ones(3), [1.0, 0.1, 2.0])
const test_settings = vec(collect(Iterators.product(cases, inits)))

const res   = map(fitter, test_settings)
serialize("fit_results.serialized", res)

