
@inline isnearint(v, tol) = abs(v-round(v)) < tol

# Unlike the previous version of this function, this inner function now ASSUMES
# that v isa Dual, and so it only hits branches that are relevant for AD.
function _besselk(v, x, maxit, tol, order)
  if abs(x) < 4 && (v < 2.85) && !isnearint(v, 0.01)
    return _besselk_ser(v, x, maxit, tol, false)
  elseif abs(x) < 8.5
    return _besselk_temme(v, x, maxit, tol, false)
  elseif abs(x) < 15.0
    return _besselk_asv(v, x, Val(12), Val(false))
  elseif abs(x) < 30.0
    return _besselk_asv(v, x, Val(8),  Val(false))
  elseif abs(v) > 1.5
    return _besselk_asv(v, x, Val(6),  Val(false))
  else
    return _besselk_as(v, x, order)
  end
end

# Just has some different cutoffs, which for whatever reason work a bit better.
# At some point this function could be improved a lot, which is part of why I'm
# okay with splitting it like this for the moment.
function _besselkxv(v, x, maxit, tol, order)
  if abs(x) < 4 && (v < 5.75) && !isnearint(v, 0.01)
    return _besselk_ser(v, x, maxit, tol, true)
  elseif abs(x) < 6.0
    return _besselk_temme(v, x, maxit, tol, true)
  elseif abs(x) < 15.0
    return _besselk_asv(v, x, Val(12), Val(true))
  elseif abs(x) < 30.0
    return _besselk_asv(v, x, Val(8), Val(true))
  elseif abs(v) > 1.5
    return _besselk_asv(v, x, Val(6), Val(true))
  else
    return _besselk_as(v, x, order)*exp(v*log(x)) # temporary, until float pows in 1.9.
  end
end

adbesselk(v::AbstractFloat, x::AbstractFloat) = Bessels.besselk(v, x)
adbesselk(v, x) = _besselk(v, x, 100, 1e-12, 6)

# TODO (cg 2022/09/09 12:22): with newer julia and/or package versions, I'm
# getting allocations in the second derivative if I'm not careful. So for now
# I'm going back to this, which unfortunately is still NaN at zero, even though
# that value is well-defined. I suppose I could put the limit in if x is zero,
# but I don't love that.
function adbesselkxv(v::AbstractFloat, x::AbstractFloat) 
  iszero(x) && return _gamma(v)*2^(v-1)
  Bessels.besselk(v, x)*(x^v)
end
adbesselkxv(v, x) = is_primal_zero(x) ? _gamma(v)*2^(v-1) : _besselkxv(v, x, 100, 1e-12, 6)

