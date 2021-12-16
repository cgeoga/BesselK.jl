module BesselK

  using LinearAlgebra, SpecialFunctions, Polynomials

  export _besselk

  include("besk_ser.jl")   # enhanced direct series. The workhorse for small-ish args.
  include("besk_as.jl")    # asymptotic expansion for large arg.
  include("besk_asv.jl")   # uniform expansion for large order.
  include("besk_temme.jl") # Temme recurrence series for small-ish args. For AD.
  include("besk.jl")       # putting it all together with appropriate branching.

end # module
