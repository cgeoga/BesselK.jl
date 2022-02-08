
# BesselK.jl

[build-latest-img]: https://github.com/cgeoga/BesselK.jl/workflows/CI/badge.svg
[build-url]: https://github.com/cgeoga/BesselK.jl/actions?query=workflow

[![][build-latest-img]][build-url] 

This package implements one function: the modified second-kind Bessel function
Kᵥ(x). It is designed specifically to be automatically differentiable **with
ForwardDiff.jl**, including providing derivatives with respect to the order
parameter `v` **that are fast and non-allocating in the entire domain for both
first and second order**.

Derivatives with respect to \nu are significantly faster than any finite
differencing method, including the most naive fixed-step minimum-order method,
and in almost all of the domain are meaningful more accurate. Particularly near
the origin you should expect to gain at least 3-5 digits. Second derivatives are
even more dramatic, both in terms of the speedup and accuracy gains, now
commonly giving 10+ more digits of accuracy.

As a happy accident/side-effect, if you're willing to give up the last couple
digits of accuracy, you could also use `ForwardDiff.jl` on this code for
derivatives with respect to argument for an order-of-magnitude speedup. In some
casual testing the argument-derivative errors with this code are never worse
than `1e-12`, and they turn 1.4 μs with allocations into 140 ns without any
allocations. 

In order to avoid naming conflicts with `SpecialFunctions.besselk`, this package
exports two functions: `adbesselk` and `adbesselkxv`. The first function is
Kᵥ(x), and the second function is (xᵛ)*Kᵥ(x). This second function has
the nice property of being bounded at the origin when v>0, and comes up in the
Matern covariance function, which was the primary motivation for this
implementation. The function `adbesselk` returns `SpecialFunctions.besselk` if
`v isa AbstractFloat`, since the AMOS `besselk` is slightly more accurate, and
there is a rule in place for the exact argument derivatives. But otherwise, it
returns `BesselK._besselk(v, x, args...)`, which is the Julia-native
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
word, this code uses several routines to evaluate Kᵥ accurately on different
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

If you use this package in your research that gets compiled into some kind of
report/article/poster/etc, please cite [this paper](https://arxiv.org/abs/2201.00090):
```
@misc{GMSS_2022,
      title={Fitting Mat\'ern Smoothness Parameters Using Automatic Differentiation}, 
      author={Christopher J. Geoga and Oana Marin and Michel Schanen and Michael L. Stein},
      year={2022},
      eprint={2201.00090},
      archivePrefix={arXiv},
      primaryClass={stat.CO}
}
```
While this package ostensibly only covers a single function, putting all of this
together and making it this fast and accurate was really a lot of work. I would
*really* appreciate you citing this paper if this package was useful in your
research. Like, for example, if you used this package to fit a Matern smoothness
parameter with second order optimization methods.

Also, if you're reading this a few months into 2022 or later, we would also
really appreciate it if you check back here or even open an issue/email to ask
if there is an official journal reference by that point. Thanks in advance!

