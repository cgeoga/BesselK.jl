
function gradient_replicates!(gstore, dfns, parms, pts, datav)
  fill!(gstore, zero(eltype(gstore)))
  @assert length(gstore) == length(dfns) "Length of gstore and number of derivative functions don't match."
  _S = Symmetric([matern(x, y, parms) for x in pts, y in pts])
  S  = cholesky(_S)
  for j in eachindex(gstore)
    dfj = dfns[j]
    Sj  = Symmetric([dfj(x, y, parms) for x in pts, y in pts])
    gstore[j] = (tr(S\Sj)*length(datav) - sum(z->dot(z, S\(Sj*(S\z))), datav))/2
  end
  gstore
end

