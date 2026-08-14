@inline isnearint(v, tol) = abs(v-round(v)) < tol

# Unlike the previous version of this function, this inner function now ASSUMES
# that v isa Dual, and so it only hits branches that are relevant for AD.
function _besselk(v, x, maxit, tol, order)
  if abs(x) < 0.5 && (v < 2.85) && !isnearint(v, 0.01)
    return _besselk_ser(v, x, maxit, tol, false)
  elseif abs(x) < 0.5
    return _besselk_temme(v, x, maxit, tol, false)
  elseif abs(x) <= 20.0
    return _besselk_intermediate(v, x, false)
  elseif abs(x) < 30.0
    return _besselk_asv(v, x, Val(8), Val(false))
  elseif abs(v) > 1.5
    return _besselk_asv(v, x, Val(6), Val(false))
  else
    return _besselk_as(v, x, order)
  end
end

# Just has some different cutoffs, which for whatever reason work a bit better.
# At some point this function could be improved a lot, which is part of why I'm
# okay with splitting it like this for the moment.
function _besselkxv(v, x, maxit, tol, order)
  if abs(x) < 0.5 && (v < 5.75) && !isnearint(v, 0.01)
    return _besselk_ser(v, x, maxit, tol, true)
  elseif abs(x) < 0.5
    return _besselk_temme(v, x, maxit, tol, true)
  elseif abs(x) <= 20.0
    return _besselk_intermediate(v, x, true)
  elseif abs(v) > 1.5
    return _besselk_asv(v, x, Val(8), Val(true))
  elseif abs(x) < 30.0
    return _besselk_asv(v, x, Val(8), Val(true))
  else
    return _besselk_as(v, x, order)*exp(v*log(x)) # temporary, until float pows in 1.9.
  end
end

@inline function _use_centered_small_x(v::Float64, x::Float64)
  av = abs(v)

  isfinite(av) &&
    isfinite(x) &&
    0.0 < x <= _BESSELK_SMALL_X_CENTERED_MAX_X &&
    av <= _BESSELK_SMALL_X_CENTERED_MAX_ORDER &&
    av != Float64(unsafe_trunc(Int, av))
end

@inline function adbesselk(v::Float64, x::Float64)
  _use_centered_small_x(v, x) ?
    _besselk_small_centered(v, x) :
    Bessels.besselk(v, x)
end
adbesselk(v::Float32, x::Float32) = Bessels.besselk(v, x)
adbesselk(v, x) = _besselk(v, x, 100, 1e-12, 6)
adbesselkxv(v, x) = is_primal_zero(x) ? _gamma(v)*2^(v-1) : _besselkxv(v, x, 100, 1e-12, 6)

