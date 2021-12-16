
using Ipopt, Printf, Serialization

# You need to run fit.jl before running this script.
include("shared.jl")
const res = deserialize("fit_results.serialized")

function printf_sym_matrix(M::Matrix{Float64})
  (s1, s2) = size(M)
  for j in 1:(s1-1)
    for k in 1:(s2-1)
      if k >= j
        @printf "%1.2f & " M[j,k]
      else
        @printf "\\cdot & "
      end
    end
    @printf "%1.2f\\\\ \n" M[j,end]
  end
  for k in 1:(s2-1)
    @printf "\\cdot & "
  end
  @printf "%1.2f" M[end,end]
end

function write_matrix_to_file(fname, M, prefix="../../../../manuscript/tables/")
  open(prefix*fname, "w") do io
    redirect_stdout(io) do
      printf_sym_matrix(M)
    end
  end
end

#=
# look at hess vs e-fish at initializer:
write_matrix_to_file("fdefish_at_ones.tex", fishfd(ones(3)))
write_matrix_to_file("fdhess_at_ones.tex",  hessfd(ones(3)))
write_matrix_to_file("adhess_at_ones.tex",  _nllh(ones(3)))
write_matrix_to_file("refhess_at_ones.tex",  high_fd_hessian(ones(3)))

# look at hess vs e-fish at MLE:
write_matrix_to_file("fdefish_at_mle.tex",  fishfd(res[3][4]))
write_matrix_to_file("fdhess_at_mle.tex",   hessfd(res[3][4]))
write_matrix_to_file("adhess_at_mle.tex",   _nllh(res[3][4]))
write_matrix_to_file("refhess_at_mle.tex",  high_fd_hessian(res[3][4]))
=#
