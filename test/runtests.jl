
using Test, BenchmarkTools, BesselK, SpecialFunctions, FiniteDifferences, ForwardDiff

const VGRID   = range(0.25, 10.0, length=100)
const XGRID   = range(0.0, 50.0,  length=201)[2:end]
const VX      = collect(Iterators.product(VGRID, XGRID))
const REF_FD1 = central_fdm(10,1)
const REF_FD2 = central_fdm(10,2)

atolfun(tru, est) = isnan(est) ? NaN : (isinf(tru) ? 0.0 : abs(tru-est))

besselkxv(v,x) = besselk(v,x)*(x^v)

fd_dbesselk_dv(v, x) = REF_FD1(_v->besselk(_v, x), v)
fd2_dbesselk_dv_dv(v, x) = REF_FD2(_v->besselk(_v, x), v)

fd_dbesselkxv_dv(v, x) = REF_FD1(_v->besselkxv(_v, x), v)
fd2_dbesselkxv_dv_dv(v, x) = REF_FD2(_v->besselkxv(_v, x), v)

ad_dbesselk_dv(v, x) = ForwardDiff.derivative(_v->BesselK.adbesselk(_v, x), v)
ad2_dbesselk_dv_dv(v, x) = ForwardDiff.derivative(_v->ad_dbesselk_dv(_v, x), v)

ad_dbesselkxv_dv(v, x) = ForwardDiff.derivative(_v->BesselK.adbesselkxv(_v, x), v)
ad2_dbesselkxv_dv_dv(v, x) = ForwardDiff.derivative(_v->ad_dbesselkxv_dv(_v, x), v)

# direct accuracy:
@testset "direct eval" begin
  println("\nDirect evaluations:")
  for (ref_fn, cand_fn, case) in ((besselk,   BesselK._besselk, :standard),
                                  (besselkxv, BesselK.adbesselkxv, :rescaled))
    amos_ref  = map(vx->ref_fn(vx[1], vx[2]), VX)
    candidate = map(vx->cand_fn(vx[1], vx[2]), VX)
    atols     = map(a_c->atolfun(a_c[1], a_c[2]), zip(amos_ref, candidate))
    ix        = findall(x-> x <= 1000.0, amos_ref)
    thresh    = case == :standard ? 5e-11 : 2e-12
    (maxerr, maxix) = findmax(abs, atols[ix])
    (maxerr_v, maxerr_x) = VX[ix][maxix]
    println("Case $case:")
    println("worst (v,x):  ($maxerr_v, $maxerr_x)")
    println("Ref value:     $(amos_ref[ix][maxix])")
    println("Est value:     $(candidate[ix][maxix])")
    println("Abs error:     $(round(maxerr, sigdigits=3))")
    @test maxerr < thresh
  end
  println()
end

# test derivative accuracy:
@testset "first derivative" begin
  println("\nFirst derivatives:")
  for (ref_fn, cand_fn, case) in ((fd_dbesselk_dv, ad_dbesselk_dv, :standard),
                                  (fd_dbesselkxv_dv, ad_dbesselkxv_dv, :rescaled))
    amos_ref  = map(vx->ref_fn(vx[1], vx[2]), VX)
    candidate = map(vx->cand_fn(vx[1], vx[2]), VX)
    atols     = map(a_c->atolfun(a_c[1], a_c[2]), zip(amos_ref, candidate))
    ix        = findall(x-> x <= 1000.0, amos_ref)
    thresh    = case == :standard ? 4e-9 : 8e-7
    (maxerr, maxix) = findmax(abs, atols[ix])
    (maxerr_v, maxerr_x) = VX[ix][maxix]
    println("Case $case:")
    println("worst (v,x):  ($maxerr_v, $maxerr_x)")
    println("Ref value:     $(amos_ref[ix][maxix])")
    println("Est value:     $(candidate[ix][maxix])")
    println("Abs error:     $(round(maxerr, sigdigits=3))")
    @test maxerr < thresh
  end
  println()
end

# test second derivative accuracy:
@testset "second derivative" begin
  println("\nSecond derivatives:")
  for (ref_fn, cand_fn, case) in ((fd2_dbesselk_dv_dv, ad2_dbesselk_dv_dv, :standard),
                                  (fd2_dbesselkxv_dv_dv, ad2_dbesselkxv_dv_dv, :rescaled))
    amos_ref  = map(vx->ref_fn(vx[1], vx[2]), VX)
    candidate = map(vx->cand_fn(vx[1], vx[2]), VX)
    atols     = map(a_c->atolfun(a_c[1], a_c[2]), zip(amos_ref, candidate))
    ix        = findall(x-> x <= 100.0, amos_ref)
    thresh    = case == :standard ? 5e-7 : 2e-6
    (maxerr, maxix) = findmax(abs, atols[ix])
    (maxerr_v, maxerr_x) = VX[ix][maxix]
    println("Case $case:")
    println("worst (v,x):  ($maxerr_v, $maxerr_x)")
    println("Ref value:     $(amos_ref[ix][maxix])")
    println("Est value:     $(candidate[ix][maxix])")
    println("Abs error:     $(round(maxerr, sigdigits=3))")
    @test maxerr < thresh
  end
  println()
end

# Testing the _xv versions really slows down the test script, and in general
# there are no no routines.
@testset "confirm no allocations" begin
  VGRID_ALLOC = (0.25, 1.0-1e-8, 1.0, 1.5, 2.1, 3.0, 3.5, 4.8)
  XGRID_ALLOC = range(0.0, 50.0, length=11)[2:end] 
  VX_ALLOC    = collect(Iterators.product(VGRID_ALLOC, XGRID_ALLOC))
  ad_alloc_test(v,x)     = @ballocated ad_dbesselk_dv($v,$x) samples=5
  ad2_alloc_test(v,x)    = @ballocated ad2_dbesselk_dv_dv($v,$x) samples=5
  ad_alloc_test_xv(v,x)  = @ballocated ad_dbesselkxv_dv($v,$x) samples=5
  ad2_alloc_test_xv(v,x) = @ballocated ad2_dbesselkxv_dv_dv($v,$x) samples=5
  ad_allocs  = map(vx->ad_alloc_test(vx[1], vx[2]), VX_ALLOC)
  ad2_allocs = map(vx->ad2_alloc_test(vx[1], vx[2]), VX_ALLOC)
  ad_allocs_xv  = map(vx->ad_alloc_test_xv(vx[1], vx[2]), VX_ALLOC)
  ad2_allocs_xv = map(vx->ad2_alloc_test_xv(vx[1], vx[2]), VX_ALLOC)
  @test all(iszero, ad_allocs)
  @test all(iszero, ad2_allocs)
  @test all(iszero, ad_allocs_xv)
  @test all(iszero, ad2_allocs_xv)
end

