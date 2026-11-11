
function pair_orders(v)
  mu = abs(v)
  mu <= 1 + 1e-3 && return (mu, mu, 0)
  n     = Int(floor(mu))
  alpha = mu - n
  (alpha, alpha, n)
end

function forward_recurrence(km, k, alpha, n::Int, x)
  for j in 0:(n-1)
    nu = alpha + j
    kn = km + 2*nu*k/x
    km = k
    k  = kn
  end
  k
end

function adbesselk(v, x)
  (isnan(x) || isnan(v) || x < zero(x)) && return NaN
  (a, alpha, n) = pair_orders(v) 
  (k_a, k_am1) = if x <= 2.0
    besselk_small_pair(a, x)
  elseif x <= 20.0
    besselk_intermediate_pair(a, x)
  else
    besselk_large_pair(a, x)
  end
  forward_recurrence(k_am1, k_a, alpha, n, x)
end

adbesselkxv(v, x) = x^v * adbesselk(v, x)

