
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

open("uk_polys.jl", "w") do out
  uk_polys = Uk_polynomials(20)
  names    = String[]
  redirect_stdout(out) do 
    run(`cat ukpoly.jl`)
    println("\n\n\n")
    for (j,pj) in enumerate(uk_polys)
      c = pj.coeffs
      stem = string("uk_", j-1)
      pnm  = string(stem, "_poly")
      push!(names, string(pnm, ","))
      str1 = string("const ", pnm, "=UkPolynomial([", map(x->string(x, ", "), c)..., "])")
      println(str1)
    end
    println()
    println(string("const UK_POLYS = [", reduce(*, names), "]"))
    println()
  end
end

