using Plots
using DelimitedFiles

const HOME = homedir()
const PROJDIR = "$HOME/Documents/Quantum_Information"
const PLOTSDIR = "$PROJDIR/graphs"
dim = 4

NAME = "exp"
datafile = "$PROJDIR/random_matrices/PS$(dim)stats_$(NAME)2.csv"
plotdata = readdlm(datafile, ',')
tolerance = 0.01
successes = plotdata[(plotdata[:,2]./plotdata[:,4]).<tolerance, :]

succ_rate = size(successes)[1]/size(plotdata)[1]

print("The success rate is $succ_rate")

fullplot = plot(plotdata[:,3], plotdata[:, 4], seriestype=:scatter,
                xlabel="NL angle", ylabel="Probability",
                legend=:none, markersize=1.5, title="all_$NAME",
                xlims=(-pi/2, pi/2), ylims=(0,1))
plot!(x->abs2(cos(2*x)+3)/16)

display(fullplot)

succplot = plot(successes[:,3], successes[:, 4], seriestype=:scatter,
                xlabel="NL angle", ylabel="Probability",
                legend=:none, markersize=1.5, title="succ_$NAME",
                xlims=(-pi/2, pi/2), ylims=(0,1))
plot!(succplot, x->abs2(cos(2*x)+3)/16)

savefig(succplot, "$PLOTSDIR/succ_$(NAME)_$dim.svg")
