using Plots
using DelimitedFiles
import JSON
 
TOLERANCE=0.01

"""Filter the data by selecting the data points for which the
optimization conververged better than TOLERANCE."""
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

function succ_filter(data, succ_data, relat_data; tol=1)
    successes = data[(succ_data./relat_data).<tol, :]
    succ_rate = size(successes)[1]/size(data)[1]
    println("The success rate is $succ_rate.")
    return successes
end

function succ_filter(data, succ_index::Integer, relat_index::Integer; tol=1)
    succ_data = view(data, :, succ_index)
    relat_data = view(data, :, relative)
    succ_filter(data, succ_data, relative=relat_data; tol=tol)
end

"""Convinient shorthand for drawaing our plot."""
function dataplot(xdata, ydata...; plot_title, width=pi/2)
    the_plot = plot!(x->abs2(cos(2*x)+3)/16) #Plot the theoretical limit curve.
    for sdata in ydata
    the_plot = plot!(xdata, sdata, seriestype=:scatter,
                     xlabel="NL angle", ylabel="Probability",
                     legend=:none, markersize=1.5, title=plot_title,
                     xlims=(-width, width), ylims=(0,1),
                     )
    end
    return the_plot
end

function quickplot(params, filepath; yindex=3, xindex=4, 
                                     width=pi/2, format = "svg")
    data = succ_filter(readdlm(filepath, ',', skipstart=1), 2)
    succplot = dataplot(data[:,xindex], data[:,yindex];
                        plot_title=params.plot_title, width=width)
    plotpath = "$(params.plot_dir)/succ_$(params.id)_$(params.plot_title).$format"
    savefig(succplot, plotpath)
    println("Plot saved to \'$plotpath\'.")
end

dicttonamedtuple(dic::AbstractDict) = (; zip(Symbol.(keys(dic)),values(dic))...)
jsonimport(paramsfile) = dicttonamedtuple.(JSON.parsefile(paramsfile))[2]

"""Copy-pasted form pss_importandrun.jl if you change it automatization
breaks."""
statspath(dir, dim, keyword, uid)="$dir/random_matrices/ps$(dim)stats_$keyword$uid.csv"

function quickplot(jsonfile; kwargs...)
    model, filedata = dicttonamedtuple.(JSON.parsefile(jsonfile))
    datapath = statspath(filedata.proj_dir,
                         model.nmodes,
                         model.merit_keyword,
                         filedata.id)
    quickplot(filedata, datapath; kwargs...)
end

if !isinteractive() && (abspath(PROGRAM_FILE) == @__FILE__)
    quickplot(ARGS[1])
end

