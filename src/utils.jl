# Return the ordinary value used for discrete branch and loop decisions.
# ForwardDiff adds the Dual method in the package extension.
@inline _primal(x::Real) = x

@inline function _cheby_clenshaw(x, coefs::NTuple{N,T}) where {N,T}
  bkp1 = zero(x*coefs[1])
  bkp2 = bkp1
  twox = x+x

  @inbounds for j in N:-1:2
    bk = muladd(twox, bkp1, coefs[j]) - bkp2
    bkp2 = bkp1
    bkp1 = bk
  end

  muladd(x, bkp1, coefs[1]) - bkp2
end

@inline function _cheby2d(sx, sy, coefs::NTuple{N,T}) where {N,T}
  inner = ntuple(i->_cheby_clenshaw(sx, coefs[i]), Val(N))
  _cheby_clenshaw(sy, inner)
end
