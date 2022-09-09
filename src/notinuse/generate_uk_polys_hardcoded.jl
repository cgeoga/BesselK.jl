
using Polynomials

# If you don't do the branch at j == 9 or so, you will totally blow up
# compilation time. Really trying to keep that under control while still
# squeezing out all allocations.
#
@inline _eta(z) = sqrt(1+z*z)+log(z/(1+sqrt(1+z*z)))
@inline _p(z)   = inv(sqrt(1+z*z))

function Uk_polynomials(max_order)
  P0   = Polynomial([1.0])
  out  = [P0]
  mul_int  = Polynomial([1.0, 0.0, -5.0])
  mul_frnt = Polynomial([0.0, 0.0, 1.0, 0.0, -1.0])/2
  for j in 1:max_order
    Pjm1 = out[end]
    Pjm1_int = integrate(mul_int*Pjm1)/8
    Pjm1_drv = derivative(Pjm1)
    newP     = mul_frnt*Pjm1_drv + Pjm1_int - Pjm1_int(0.0)/8
    push!(out, newP)
  end
  out
end

open("uk_polys_hardcoded.jl", "w") do out
  uk_polys = Uk_polynomials(20)
  redirect_stdout(out) do 
    for (j,pj) in enumerate(uk_polys)
      c = pj.coeffs
      if j <= 10
        str = string("uk_", j-1, "_eval(x)= @evalpoly x", map(x->string(" ", x, " "), c)...)
        println(str)
      else
        stem = string("uk_", j-1)
        str1 = string("const ", stem, "_coef=(", map(x->string(x, ", "), c)..., ")")
        println(str1)
        cnm  = string(stem, "_coef")
        str2 = string(stem, "_eval(x) = evalpoly(x, $cnm)")
        println(str2)
      end
      println()
    end
    str = string(
    "const UK_POLYS = (uk_0_eval, uk_1_eval, uk_2_eval, uk_3_eval, uk_4_eval,
    uk_5_eval, uk_6_eval, uk_7_eval, uk_8_eval, uk_9_eval, uk_10_eval,
    uk_11_eval, uk_12_eval, uk_13_eval, uk_14_eval, uk_15_eval, uk_16_eval,
    uk_17_eval, uk_18_eval, uk_19_eval, uk_20_eval)")
    println(str)
  end
end

