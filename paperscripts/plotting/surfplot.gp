
set terminal epslatex color blacktext
#set terminal pngcairo size 800,500 enhanced font 'Verdana,15'

set key off

load 'cmocean_deep.pal'

set datafile separator ","
file="../plotdata/profile_surface.csv"

set cbrange[0:12]

set autoscale xfix
set autoscale yfix

set xlabel '$\nu$'
set ylabel '$\rho$'

set output "profile_liksurface.tex"

plot file matrix nonuniform with image, \
     "contours.dat" using 1:2 "%lf %lf %lf" with lines lc rgb "black"

unset output

