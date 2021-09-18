using Plots
using DelimitedFiles

const HOME = homedir()
PARAMETERS_FILE="pss_parameter_files/dim3-ph2-a.jl" #File containing execution parameters
include(PARAMETERS_FILE)

"Can this be rewritten to support other kinds of index?"
function succ_filter(data, succ_data; tol=TOLERANCE)
    successes = data[succ_data.<tol, :]
    succ_rate = size(successes)[1]/size(data)[1]
    print("The success rate is $succ_rate")
    return successes
end


function succ_filter(data, succ_index::Integer; tol=TOLERANCE)
    succ_data = view(data, :, succ_index)
    succ_filter(data, succ_data)
end

function succ_filter(data, succ_data, relat_data; tol=TOLERANCE)
    successes = data[(succ_data./relat_data).<tol, :]
    succ_rate = size(successes)[1]/size(data)[1]
    print("The success rate is $succ_rate.\n")
    return successes
end

function succ_filter(data, succ_index::Integer, relat_index::Integer; tol=TOLERANCE)
    succ_data = view(data, :, succ_index)
    relat_data = view(data, :, relative)
    succ_filter(data, succ_data, relative=relat_data)
end

function dataplot(xdata, ydata...; plot_title=PLOT_TITLE) 
    the_plot = plot!(x->abs2(cos(2*x)+3)/16) #Plot the theoretical limit curve.
    for sdata in ydata
    the_plot = plot!(xdata, sdata, seriestype=:scatter,
                     xlabel="NL angle", ylabel="Probability",
                     legend=:none, markersize=1.5, title=plot_title,
                     xlims=(-pi/2, pi/2), ylims=(0,1),
                     )
    end
    return the_plot           
end

function quickplot(ang_index, prob_index, datafile=IMPNAME)
    data = succ_filter(readdlm(datafile, ',', skipstart=1), 2) 
    succplot = dataplot(data[:,ang_index], data[:,prob_index])
    savefig(succplot, "$PLOTSDIR/succ_$(NAME)_$DIM.svg")
    print("Plot saved to $PLOTSDIR/succ_$(NAME)_$DIM.svg\n")
end

if !isinteractive()
    quickplot(4,3)
end