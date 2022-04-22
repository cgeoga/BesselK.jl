mkdir -p ../plotdata/
mkdir -p ../plotdata/tables/

$JULIA_HOME/julia --project accuracy.jl
$JULIA_HOME/julia --project derivative_accuracy.jl
$JULIA_HOME/julia --project secderivative_accuracy.jl
$JULIA_HOME/julia --project speed1.jl >speed_results.tex
$JULIA_HOME/julia --project speed2.jl >matern_results.tex

mv speed_results.tex  ../plotdata/tables/
mv matern_results.tex ../plotdata/tables/
