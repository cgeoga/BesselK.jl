
module BesselK

  using Bessels

  export adbesselk, adbesselkxv, matern

  # Here's a work-around since ForwardDiff#481 got merged, making iszero(x)
  # check that the value AND partials of x are zero. Conceptually, I'm
  # sympathetic that this is the more correct choice. It just doesn't quite work
  # for the way this code needs to branch.
  _iszero(x) = ifelse(0.0 <= x <= 0.0, true, false)

  include("gamma.jl")      # gamma function, for the moment ripped from Bessels.jl
  include("besk_ser.jl")   # enhanced direct series. The workhorse for small-ish args.
  include("besk_as.jl")    # asymptotic expansion for large arg.
  include("uk_polys.jl")   # Uk polynomials, now generated statically ahead of time.
  include("besk_asv.jl")   # uniform expansion for large order.
  include("besk_temme.jl") # Temme recurrence series for small-ish args. For AD.
  include("besk.jl")       # putting it all together with appropriate branching.
  include("matern.jl")     # a basic Matern covariance function

end 

