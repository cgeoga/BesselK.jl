
@generated function _besselk_vz_asv(v, z, ::Val{N}, ::Val{M}) where{N,M}
  quote
    ez  = sqrt(1+z*z)+log(z/(1+sqrt(1+z*z)))
    pz  = inv(sqrt(1+z*z))
    out = sqrt(pi/(2*v))/sqrt(sqrt(1+z*z))
    mulval = M ? exp(v*log(z*v)-v*ez) : exp(-v*ez)
    (ser, sgn, _v) = (zero(z), one(z), one(v))
    evaled_polys = tuple($([:($(Symbol(:uk_, j, :_poly))(pz)) for j in 0:(N-1)]...))
    Base.Cartesian.@nexprs $N j -> begin
      ser += sgn*evaled_polys[j]/_v
      sgn *= -one(z)
      _v  *= v
    end
    mulval*out*ser
  end
end

_besselk_asv(v, z, maxorder, modify) = _besselk_vz_asv(v, z/v, maxorder, modify)

