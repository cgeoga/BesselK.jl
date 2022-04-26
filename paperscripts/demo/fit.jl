
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
  is_bfgs = in(case, (:FD_BFGS, :AD_BFGS))
  println("Case: $case")
  cache = Vector{Float64}[]
  box_l = is_bfgs ? [1e-2, 1e-2, 0.25] : [0.0, 0.0, 0.25]
  prob  = createProblem(3, box_l, fill(1e22, (3,)), 0,
                        Float64[], Float64[], 0, div(3*4, 2),
                        _p->caching_nll(_p, cache),
                        (args...)->nothing,
                        gradfun,
                        (args...)->nothing,
                        (x,m,r,c,o,l,v)->ipopt_hessian(x,m,r,c,o,l,v,
                                                       hesfun,Function[],0))
  addOption(prob, "tol", 1e-5)
  addOption(prob, "max_iter", maxiter)
  if is_bfgs
    addOption(prob, "hessian_approximation", "limited-memory")
    addOption(prob, "nlp_scaling_method", "none")
  end
  prob.x = deepcopy(ini) # for safety to avoid weird persistent pointer games.
  try
    @time status = solveProblem(prob)
    _h = is_bfgs ? fil(NaN, 3, 3) : hesfun(prob.x)
    return (prob, status, case, deepcopy(prob.x), _nll(prob.x), _h, cache, ini)
  catch er
    println("\n\nOptimization failed with error $er\n\n")
    return (prob, :FAIL_OPT_ERR, case, deepcopy(prob.x), NaN, fill(NaN, 3, 3), cache, ini)
  end
end

const cases = ((:FD_FISH, grad_fd!, fishfd, 100),
               (:AD_FISH, grad_ad!, fishad, 100),
               (:FD_HESS, grad_fd!, hessfd, 100), 
               (:AD_HESS, grad_ad!, _nllh,  100),
               (:FD_BFGS, grad_fd!, no_hessian, 100), 
               (:AD_BFGS, grad_ad!, no_hessian, 100))

const inits = (ones(3), [1.0, 0.1, 2.0])
const test_settings = vec(collect(Iterators.product(cases, inits)))

if !isinteractive()
  const res   = map(fitter, test_settings)
  serialize("fit_results.serialized", res)
end

