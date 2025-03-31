
module BesselK

  using Bessels

  export adbesselk, adbesselkxv, matern

  # A special method that has an extension for ForwardDiff to handle x::Dual.
  is_primal_zero(x::Real) = iszero(x)

  include("gamma.jl")      # gamma function, for the moment ripped from Bessels.jl
  include("besk_ser.jl")   # enhanced direct series. The workhorse for small-ish args.
  include("besk_as.jl")    # asymptotic expansion for large arg.
  include("uk_polys.jl")   # Uk polynomials, now generated statically ahead of time.
  include("besk_asv.jl")   # uniform expansion for large order.
  include("besk_temme.jl") # Temme recurrence series for small-ish args. For AD.
  include("besk.jl")       # putting it all together with appropriate branching.
  include("matern.jl")     # a basic Matern covariance function

end 

