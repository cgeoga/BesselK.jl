
module BesselKForwardDiffExt

  using BesselK, ForwardDiff

  import ForwardDiff: Dual

  BesselK.is_primal_zero(x::Dual) = BesselK.is_primal_zero(ForwardDiff.value(x))
  BesselK._primal(x::Dual) = BesselK._primal(ForwardDiff.value(x))

end

