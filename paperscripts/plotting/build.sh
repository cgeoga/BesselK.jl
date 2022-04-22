
mkdir -p ../plotdata/figures/

gnuplot *.gp
mv *.eps ../plotdata/figures/
mv *.tex ../plotdata/figures/

