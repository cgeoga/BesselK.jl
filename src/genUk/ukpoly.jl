
struct UkPolynomial{N,P}
  coef_skip_zeros::P
  constant::Float64
end

function UkPolynomial(coefv::Vector{T}) where{T}
  # find the first non-zero coefficient, assuming the first coefficient is a
  # zero-order constant:
  first_nonzero_index = findfirst(!iszero, coefv)
  # now take the coefficients that are not zero, which is every other one:
  if first_nonzero_index == 1
    constant = coefv[1]
    next_nonzero_index  = findfirst(!iszero, coefv[2:end])
    if isnothing(next_nonzero_index)
      if length(coefv) > 1 
        throw(error("These coefficients don't correspond to a U_k polynomial."))
      end
      return UkPolynomial{0,Nothing}(nothing, coefv[1])
    end
    first_nonzero_index = next_nonzero_index+1
  else
    constant = 0.0 
  end
  nzcoef = vcat(zero(T), coefv[first_nonzero_index:2:end])
  coef_skip_zeros = length(nzcoef) < 20 ? tuple(nzcoef...) : nzcoef
  (N,P) = (first_nonzero_index-1, typeof(coef_skip_zeros))
  UkPolynomial{N,P}(coef_skip_zeros, constant)
end

# the case of a simple polynomial with no leading zeros.
(Uk::UkPolynomial{0,P})(x) where{P} = evalpoly(x^2, Uk.coef_skip_zeros)/x + Uk.constant

# the case of a zero-order polynomial:
(Uk::UkPolynomial{0,Nothing})(x) = Uk.constant

# the nontrivial case of a polynomial that DOES have leading zeros:
function (Uk::UkPolynomial{N,P})(x) where{N,P}
  pre_multiply = x^(N-2)
  pv = evalpoly(x^2, Uk.coef_skip_zeros)
  pre_multiply*pv + Uk.constant
end

