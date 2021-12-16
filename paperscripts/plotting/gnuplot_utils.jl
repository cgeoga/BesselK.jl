
function tablesave!(name, M::Matrix; rlabels=nothing, clabels=nothing, siunits=false)
  out = siunits ? map(x->"\\num{"*string(x)*"}", M) : string.(M)
  if !isnothing(rlabels)
    out = hcat(string.(rlabels), out)
  end
  if !isnothing(clabels)
    is_row_nothing = isnothing(rlabels)
    rw1 = string.(clabels)
    if !is_row_nothing
      pushfirst!(rw1, "")
    end
    rw1[end] = rw1[end] * "\\\\"
    writedlm("header_"*name, reshape(rw1, 1, length(rw1)), '&')
  end
  out[:,end] .= out[:,end] .* repeat(["\\\\"], size(out,1))
  writedlm(name, out, '&')
end

function gnuplot_save_matrix!(name, M::Matrix{Float64}, row_pts, col_pts, delim=',')
  out = Array{Any}(undef, size(M,1)+1, size(M,2)+1)
  out[1,1] = ""
  out[1,2:end] .= col_pts
  out[2:end,1] .= row_pts
  out[2:end, 2:end] .= M
  writedlm(name, out, delim)
end

function gnuplot_save_vector!(name, M, pts, delim=',')
  writedlm(name, hcat(pts, M), delim)
end

