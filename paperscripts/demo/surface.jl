
include("shared.jl")
include("../plotting/gnuplot_utils.jl")

BLAS.set_num_threads(1)

# Small neighborhood of the MLE.
const RANGE_GRID  = range(1.5,  6.0,  length=60)
const SMOOTH_GRID = range(1.2, 1.4, length=50)

function gensurf()
  out = zeros(length(RANGE_GRID), length(SMOOTH_GRID))
  Threads.@threads for j in eachindex(RANGE_GRID)
    @inbounds for k in eachindex(SMOOTH_GRID)
      try
        out[j,k] = profile_nll_replicates((RANGE_GRID[j], SMOOTH_GRID[k]), PTS, SIMS)
      catch
        println("Warning, caught a failure.")
        out[j,k] = NaN 
      end
    end
  end
  out .- minimum(out)
end

# Will need to figure out the right color scale here.
const surf = gensurf()

gnuplot_save_matrix!("../plotdata/profile_surface.csv", surf, 
                     RANGE_GRID, SMOOTH_GRID)

