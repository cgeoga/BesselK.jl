
# BesselK.jl

**NOTE:** As of release `v0.5.9`, this library is a Julia port of
[`matern`](https://git.sr.ht/~cgeoga/libmatern-rs), a Rust crate that I
developed using new routines based on the old source code of this library and of
[`Bessels.jl`](https://github.com/JuliaMath/Bessels.jl).

[build-latest-img]: https://github.com/cgeoga/BesselK.jl/workflows/CI/badge.svg
[build-url]: https://github.com/cgeoga/BesselK.jl/actions?query=workflow

[![][build-latest-img]][build-url] 

This package implements the following functions:

- `adbesselk(v, x)`: the modified second-kind Bessel function `besselk`, except
  with strategies that are entirely auto-differentiable with respect to order as
  well as argument.

- `adbesselkxv(v, x)`: for now, a very simple `x^v * adbesselk(v, x)`. But in
  the future, more specialized routines may be developed.

- `matern(dist, params)`: the [Matern covariance function](https://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function)
  parameterized with `(\sigma, \rho, \nu)`. There is also a `matern(x, y, params)`
  that evaluates to `matern(norm(x-y), params)`.

# Limitations

* **Only derivatives up to the second are checked and confirmed accurate.** The
  code uses a large number of local polynomial expansions at slightly hairy
  values of internal intermediate functions, and so at some sufficiently high
  level of derivative those local polynomials won't give accurate partial
  information.

# Citation

If you use this package in your research that gets compiled into some kind of
report/article/poster/etc, please cite [this paper](https://arxiv.org/abs/2201.00090):
```
@misc{GMSS_2022,
      title={Fitting Mat\'ern Smoothness Parameters Using Automatic Differentiation}, 
      author={Christopher J. Geoga and Oana Marin and Michel Schanen and Michael L. Stein},
      year={2022},
      journal={Statistics and Computing}
}
```
While this package ostensibly only covers a single function, putting all of this
together and making it this fast and accurate was really a lot of work. I would
*really* appreciate you citing this paper if this package was useful in your
research. Like, for example, if you used this package to fit a Matern smoothness
parameter with second order optimization methods.

