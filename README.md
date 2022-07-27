
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

# A note to people coming here from the paper
You'll see that this repo defines a great deal of specific derivative functions
in the files in `./examples` and `./paperscripts`. **This is only because we
specifically tested those quantities in the paper**. If you're just here to fit
a Matern covariance function, then you should **not** be doing that. Your code,
at least in the simplest case, should probably look more like this:
```julia

using ForwardDiff, BesselK

function my_covariance_function(loc1, loc2, params)
  ... # your awesome covariance function, presumably using adbesselk somewhere.
end

const my_data = ... # load in your data
const my_locations = ... # load in your locations

# Create your likelihood and use ForwardDiff for the grad and Hessian:
function nll(params)
  K = cholesky!(Symmetric([my_covariance_function(x, y, params) 
                           for x in my_locations, y in my_locations]))
  0.5*(logdet(K) + dot(my_data, K\my_data))
end
nllg(params) = ForwardDiff.gradient(nll, params)
nllh(params) = ForwardDiff.hessian(nll, params)

my_mle = some_optimizer(init_params, nll, nllg, nllh, ...)
```
Or something like that. You of course do not *have* to do it this way, and could
manually implement the gradient and Hessian of the likelihood after manually
creating derivatives of the covariance function itself (see
`./example/matern.jl` for a demo of that), and manual implementations,
particularly for the Hessian, will be faster if they are thoughtful enough. But
what I mean to emphasize here is that in general you should *not* be doing
manual chain rule or derivative computations of your covariance function itself.
Let the AD handle that for you and enjoy the power that Julia's composability
offers.

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

# Also consider: `Bessels.jl`

This software package was written with the pretty specific goal of computing
derivatives of Kᵥ(x) with respect to the order using `ForwardDiff.jl`. While it
is in general a bit faster than AMOS, we give up a few digits of accuracy here
and there in the interest of better and faster derivatives. If you just want the
fastest possible Kᵥ(x), then you would probably be better off using
[`Bessels.jl`](https://github.com/heltonmc/Bessels.jl). At the time of writing
it only offers Kᵥ(x) for integer orders, but non-integer orders will be
available soon enough I'm sure. While differentiability is on the roadmap,
they more explicitly target writing the fastest possible base Kᵥ(x), and what
they offer is seriously fast.

There is and will be some cross-pollination between the two software projects,
and at some point I expect to switch `adbesselk` to use `Bessels.besselk` where
possible instead of `SpecialFunctions.besselk`. And  at some point if the order
derivatives become available there might not be much reason to use this package
instead of that one, although I think for the moment if you want to fit a Matern
covariance function you probably need to be here.

On the topic, the following methods are lifted directly from `Bessels.jl` so
that we can go fast in the meantime:

* Integer order `\nu` when `\nu isa AbstractFloat`.

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

A to-do item (written 2022/07/27), I think, is to re-organize the code a bit so
that there is a function `_besselk_vdual` that only gets called when `v isa
ForwardDiff.Dual` and a function `_besselk_abstractfloat` when `v isa
AbstractFloat`. For the initial release, `adbesselk` always defaulted to
`SpecialFunctions.besselk` where possible to give people what they expected and
every digit possible. But as `Bessels.jl` matures, I think lifting at least a
few of those routines in the interim is appealing but means that there is
awkwardly a lot of control flow in `BesselK._besselk` as well as
`BesselK.adbesselk` now. Probably better to compartmentalize those two domain
partitionings. Defaulting to `Bessels.besselk` when `v isa
AbstractFloat` is probably a good intermediate goal once it's ready for all
arguments.

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

