
# Conclusion: for e-fish, pretty similar.

function efish_replicates(dfns, parms, pts, lendatav)
  fish = zeros(3,3)
  _S   = Symmetric([matern(x, y, parms) for x in pts, y in pts])
  S    = cholesky!(_S)
  for j in 1:3
    dfj = dfns[j]
    Sj  = Symmetric([dfj(x, y, parms) for x in pts, y in pts])
    fish[j,j] = tr(S\(Sj*(S\Sj)))/2
    for k in (j+1):3
      dfk = dfns[k]
      Sk  = Symmetric([dfk(x, y, parms) for x in pts, y in pts])
      fish[j,k] = tr(S\(Sj*(S\Sk)))/2
      fish[k,j] = fish[j,k]
    end
  end
  fish*lendatav
end

