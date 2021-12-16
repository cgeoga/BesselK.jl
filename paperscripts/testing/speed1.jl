
using BesselK, BenchmarkTools, Printf

include("shared.jl")

const PAIRS = ((0.5,   1.0,  "half-integer"),
               (1.0,   1.0,  "whole integer"),
               (3.001, 1.0,  "near-integer order"),
               (3.001, 8.0,  "near-integer order, borderline arg"),
               (1.85,  1.0,  "small order"),
               (1.85,  8.0,  "small order, borderline arg"),
               (1.85,  14.0, "intermediate arg"),
               (1.85,  29.0, "large intermediate arg"),
               (1.85,  35.0, "large argument"))

for (j, (v, x, descriptor)) in enumerate(PAIRS)
  t_us    = (@belapsed BesselK._besselk($v, $x)   samples=1_000)*1e9 # our code
  t_am    = (@belapsed BesselK.besselk($v, $x)    samples=1_000)*1e9 # AMOS
  t_us_d  = (@belapsed adbesselkdv($v, $x)        samples=1_000)*1e9 # our code
  t_am_d  = (@belapsed fdbesselkdv($v, $x)        samples=1_000)*1e9 # AMOS
  t_us_d2 = (@belapsed adbesselkdvdv($v, $x)      samples=1_000)*1e9 # our code
  t_am_d2 = (@belapsed fdbesselkdvdv($v, $x)      samples=1_000)*1e9 # AMOS
  if j != length(PAIRS)
    @printf "(%1.3f, %1.0f) & %1.0f & %1.0f & %1.0f & %1.0f & %1.0f & %1.0f & %s \\\\\n" v x t_us t_am t_us_d t_am_d t_us_d2 t_am_d2 descriptor
  else
    @printf "(%1.3f, %1.0f) & %1.0f & %1.0f & %1.0f & %1.0f & %1.0f & %1.0f & %s\n" v x t_us t_am t_us_d t_am_d t_us_d2 t_am_d2 descriptor
  end
end


