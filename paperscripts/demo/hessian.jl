
function hessian_replicates(dfns, d2fns, parms, pts, datav)
  hess = zeros(3,3)
  _S   = Symmetric([matern(x, y, parms) for x in pts, y in pts])
  S    = cholesky(_S)
  for j in 1:3
    dfj = dfns[j]
    Sj  = Symmetric([dfj(x, y, parms) for x in pts, y in pts])
    for k in j:3
      # first derivative matrix:
      dfk   = dfns[k]
      Sk    = Symmetric([dfk(x, y, parms)  for x in pts, y in pts])
      # second derivative matrix:
      dfjk  = d2fns[j][k-j+1]
      Sjk   = Symmetric([dfjk(x, y, parms) for x in pts, y in pts])
      # Compute the complicated derivative of the qform (not efficient or
      # thoughtful, so don't use this code for real somewhere):
      dqf   = -(S\(Sk*(S\(Sj*inv(_S)))) )
      dqf  += S\(Sjk*inv(_S))
      dqf  -= S\(Sj*(S\(Sk*inv(_S)))) 
      # compute the term:
      term  = (tr(S\Sjk) - tr(S\(Sj*(S\Sk))))*length(datav)/2 
      term -= sum(z->dot(z, dqf, z), datav)/2
      hess[j,k] = term
      hess[k,j] = term
    end
  end
  hess
end

