
using StaticArrays, BesselK, ForwardDiff

# As simple as this:
const x_point = @SVector rand(2)
const y_point = @SVector rand(2) 
const params  = [1.0, 1.0, 1.0] # scale, range, smoothness.

# All derivatives, including the smoothness!
ForwardDiff.gradient(p->matern(x_point, y_point, p), params)


