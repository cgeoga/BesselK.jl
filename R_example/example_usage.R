
# Note that you'll need to install Julia as well as some dependencies. Check out
# the github page for JuliaCall (https://github.com/Non-Contradiction/JuliaCall)
# for some documentation on how best to do that. If you really have trouble,
# please open an issue here and we can talk about what's going wrong and
# ultimately update these example files so other people don't have to work out
# the same issues from scratch.

library(JuliaCall)

julia_setup()
julia_source("R_interface.jl")

# NOTE THAT THIS TAKES (v,x), NOT THE (x,v) USED BY R. You are of course welcome
# to change this in your own adoption of this code, but then you'll have to be
# careful not to get mixed up when you go to the Julia code, which uses (v,x).
besselk_dv_dv <- function(v, x) {
  julia_call("R_besselk_dv_dv", 1.1, 2.1)
}

example_bessel = besselk_dv_dv(1.1, 1.1)
example_matern = julia_call("matern_d1_d3", 1.1, c(1.1, 1.1, 1.1))

