
using Test, BesselK, SpecialFunctions, FiniteDifferences, ForwardDiff

const VGRID   = range(0.25, 10.0, length=100)
const XGRID   = range(0.0, 50.0,  length=201)[2:end]
const VX      = collect(Iterators.product(VGRID, XGRID))
const REF_FD1 = central_fdm(10,1)
const REF_FD2 = central_fdm(10,2)

atolfun(tru, est) = isnan(est) ? NaN : (isinf(tru) ? 0.0 : abs(tru-est))

fd_dbesselk_dv(v, x) = REF_FD1(_v->besselk(_v, x), v)
ad_dbesselk_dv(v, x) = ForwardDiff.derivative(_v->BesselK.adbesselk(_v, x), v)

fd2_dbesselk_dv_dv(v, x) = REF_FD2(_v->besselk(_v, x), v)
ad2_dbesselk_dv_dv(v, x) = ForwardDiff.derivative(_v->ad_dbesselk_dv(_v, x), v)

# direct accuracy:
@testset "direct eval" begin
  amos_ref  = map(vx->besselk(vx[1], vx[2]), VX)
  candidate = map(vx->BesselK._besselk(vx[1], vx[2]), VX)
  atols     = map(a_c->atolfun(a_c[1], a_c[2]), zip(amos_ref, candidate))
  atols_f   = atols[amos_ref .<= 1000.0]
  @test maximum(atols_f) < 1e-10
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

