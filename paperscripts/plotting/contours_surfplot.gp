
set datafile separator ","

file="../plotdata/profile_surface.csv"

set pm3d map

unset surface
set contour base
set cntrparam levels incremental 0.5,1,8.5
set table "contours.dat"
splot file matrix nonuniform with lines
unset table

