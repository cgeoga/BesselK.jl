module BesselK

export adbesselk, adbesselkxv, matern

include("gamma.jl")      # gamma function, for the moment ripped from Bessels.jl
include("small_arg.jl")  # small argument routine
include("intermediate_arg.jl")  # intermediate argument routine
include("large_arg.jl")  # large argument routine
include("besk.jl")       # combined routine
include("matern.jl")     # a basic Matern covariance function

end
