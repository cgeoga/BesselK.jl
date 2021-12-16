
# BesselK.jl

This package implements one function: the modified second-kind Bessel function
K_\nu(x). It is designed specifically to be automatically differentiable **with
ForwardDiff.jl**, including providing derivatives with respect to the order
parameter \nu.

In order to avoid naming conflicts with `SpecialFunctions.besselk`, this package
exports two functions: `adbesselk` and `adbesselkxv`. The first function is
K_\nu(x), and the second function is (x^\nu)\*K_\nu(x). This second function has
the nice property of being bounded at the origin when \nu>0, and comes up in the
Matern covariance function, which was the primary motivation for this
implementation. The function `adbesselk` returns `SpecialFunctions.besselk` if
`v isa AbstractFloat`, since the AMOS `besselk` is slightly more accurate. But
otherwise, it returns `BesselK._besselk`, which is the Julia-native
implementation here that provides very accurate derivatives.

Here is a very basic demo:
```julia
using ForwardDiff, SpecialFunctions, BesselK

(v, x) = (1.1, 2.1)

# For regular evaluations, you get what you're used to getting:
@assert isapprox(besselk(v, x), adbesselk(v, x))
@assert isapprox((x^v)*besselk(v, x), adbesselkxv(v, x))

# But now you also get good (and fast!) derivatves:
@show ForwardDiff.derivative(_v->adbesselk(_v, x), v)   # good to go.
@show ForwardDiff.derivative(_v->adbesselkxv(_v, x), v) # good to go.
```

# Limitations

For the moment there are two primary limitations:

* **AD compatibility with `ForwardDiff.jl` only**. The issue here is that in one
  particular case I use a different function branch of one is taking a
  derivative with respect to `v` or just evaluating `besselk(v, x)`. The way that
  is currently checked in the code is with `if (v isa AbstractFloat)`, which may
  not work properly for other methods.

* **Only derivatives up to the second are checked and confirmed accurate.** The
  code uses a large number of local polynomial expansions at slightly hairy
  values of internal intermediate functions, and so at some sufficiently high
  level of derivative those local polynomials won't give accurate partial
  information.

# Implementation details

See the reference for an entire paper discussing the implementation. But in a
word, this code uses several routines to evaluate K_\nu accurately on different
parts of the domain, and has to use some non-standard to maintain AD
compatibility and correctness. When `v` is an integer or half-integer, for
example, a lot of additional work is required.

The code is also pretty well-optimized, and you can benchmark for yourself or
look at the paper to see that in several cases the `ForwardDiff.jl`-generated
derivatives are faster than a single call to `SpecialFunctions.besselk`. To
achieve this performance, particularly for second derivatives, some work was
required to make sure that all of the function calls are non-allocating, which
means switching from raw `Tuple`s to `Polynomial` types in places where the
polynomials are large enough and things like that. Again this arguably makes the
code look a bit disorganized or inconsistent, but to my knowledge it is all
necessary. If somebody looking at the source finds a simplification, I would
love to see it, either in terms of an issue or a PR or an email or a patch file
or anything. 

# Citation

Coming very shortly. While this package ostensibly only covers a single
function, putting all of this together and making it this fast and accurate was
really a lot of work. I would *really* appreciate you citing this paper if this
package was useful in your research. Like, for example, if you used this package
to fit a Matern smoothness parameter with second order optimization methods.

