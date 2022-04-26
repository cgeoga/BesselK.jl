
# For some reason THESE proportions make a square looking figure....
set terminal epslatex color blacktext
set key off
set datafile separator ","

load 'cmocean_balance.pal'

set autoscale xfix
set autoscale yfix
set logscale cb

#set cbtics format '\tiny %g'
#set cblabel '\footnotesize Absolute error' 

#set format x '\tiny %g'
set xlabel '$x$'

#set format y '\tiny %g'
set ylabel '$\nu$'


# Function absolute tolerances:
load 'cmocean_deep.pal'
set palette negative
set cbrange[1e-20:1e-1]
do for [case in "amos ours"] {
  set output "atols_".case.".tex"
  file="../plotdata/atols_".case.".csv"
  plot file matrix nonuniform with image
  unset output
}

# Derivative absolute tolerances:
set cbrange [1e-20:100000]
do for [case in "fd ad fdad"] {
  if (case eq "fdad"){
    load 'cmocean_balance.pal'
    set palette positive
    unset logscale cb
    set cbrange [-10:10]
  }else{
    load 'cmocean_deep.pal'
    set palette negative
    set logscale cb
  }
  set output "atols_deriv_".case.".tex"
  file="../plotdata/atols_deriv_".case.".csv"
  plot file matrix nonuniform with image
  unset output
}

# second deriv absolute tolerance:
set cbrange [1e-20:2e12]
do for [case in "fd ad fdad"] {
  if (case eq "fdad"){
    load 'cmocean_balance.pal'
    set palette positive
    unset logscale cb
    set cbrange [-10:10]
  }else{
    load 'cmocean_deep.pal'
    set palette negative
    set logscale cb
  }
  set output "atols_deriv2_".case.".tex"
  file="../plotdata/atols_deriv2_".case.".csv"
  plot file matrix nonuniform with image
  unset output
}

