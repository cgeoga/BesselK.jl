
using Test, BenchmarkTools, BesselK, SpecialFunctions, FiniteDifferences, ForwardDiff

const VGRID   = range(0.25, 10.0, length=100)
const XGRID   = range(0.0, 50.0,  length=201)[2:end]
const VX      = collect(Iterators.product(VGRID, XGRID))
const REF_FD1 = central_fdm(10,1)
const REF_FD2 = central_fdm(10,2)

atolfun(tru, est) = isnan(est) ? NaN : (isinf(tru) ? 0.0 : abs(tru-est))

besselkxv(v,x) = besselk(v,x)*(x^v)

fd_dbesselk_dv(v, x) = REF_FD1(_v->besselk(_v, x), v)
ad_dbesselk_dv(v, x) = ForwardDiff.derivative(_v->BesselK.adbesselk(_v, x), v)

fd2_dbesselk_dv_dv(v, x) = REF_FD2(_v->besselk(_v, x), v)
ad2_dbesselk_dv_dv(v, x) = ForwardDiff.derivative(_v->ad_dbesselk_dv(_v, x), v)

# direct accuracy:
@testset "direct eval" begin
  for (ref_fn, cand_fn) in ((besselk,   BesselK._besselk),
                            (besselkxv, BesselK.adbesselkxv))
    amos_ref  = map(vx->ref_fn(vx[1], vx[2]), VX)
    candidate = map(vx->cand_fn(vx[1], vx[2]), VX)
    atols     = map(a_c->atolfun(a_c[1], a_c[2]), zip(amos_ref, candidate))
    atols_f   = atols[amos_ref .<= 1000.0]
    @test maximum(atols_f) < 1e-10
  end
end

# test derivative accuracy:
@testset "first derivative" begin
  amos_ref  = map(vx->fd_dbesselk_dv(vx[1], vx[2]), VX)
  candidate = map(vx->ad_dbesselk_dv(vx[1], vx[2]), VX)
  atols     = map(a_c->atolfun(a_c[1], a_c[2]), zip(amos_ref, candidate))
  atols_f   = atols[amos_ref .<= 1000.0]
  @test maximum(atols_f) < 1e-8
end

# test second derivative accuracy:
@testset "second derivative" begin
  amos_ref  = map(vx->fd2_dbesselk_dv_dv(vx[1], vx[2]), VX)
  candidate = map(vx->ad2_dbesselk_dv_dv(vx[1], vx[2]), VX)
  atols     = map(a_c->atolfun(a_c[1], a_c[2]), zip(amos_ref, candidate))
  atols_f   = atols[amos_ref .<= 100.0]
  @test maximum(atols_f) < 1e-6
end

@testset "confirm no allocations" begin
  VGRID_ALLOC = (0.25, 1.0-1e-8, 1.0, 1.5, 2.1, 3.0, 3.5, 4.8)
  XGRID_ALLOC = range(0.0, 50.0, length=11)[2:end] 
  VX_ALLOC    = collect(Iterators.product(VGRID_ALLOC, XGRID_ALLOC))
  ad_alloc_test(v,x)  = @ballocated ad_dbesselk_dv($v,$x) samples=10
  ad2_alloc_test(v,x) = @ballocated ad2_dbesselk_dv_dv($v,$x) samples=10
  ad_allocs  = map(vx->ad_alloc_test(vx[1], vx[2]), VX_ALLOC)
  ad2_allocs = map(vx->ad2_alloc_test(vx[1], vx[2]), VX_ALLOC)
  @test all(iszero, ad_allocs)
  @test all(iszero, ad2_allocs)
end

