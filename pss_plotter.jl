using Plots
using DelimitedFiles
import JSON

TOLERANCE=0.01

"""
Filter the data by selecting the data points for which the
optimization conververged to a unitary post-selection process better
than the tolerance 'tol' .
"""
function
succ_filter(data, succ_data; tol=TOLERANCE)
    successes = data[succ_data.<tol, :]
    succ_rate = size(successes)[1]/size(data)[1]
    print("The success rate is $succ_rate")
    return successes
end

function succ_filter(data, succ_index::Integer; tol=TOLERANCE)
    succ_data = view(data, :, succ_index)
    succ_filter(data, succ_data)
end

function succ_filter(data, succ_data, relat_data; tol)
    successes = data[(succ_data./relat_data).<tol, :]
    succ_rate = size(successes)[1]/size(data)[1]
    println("The success rate is $succ_rate.")
    return successes
end

function succ_filter(data, succ_index::Integer, relat_index::Integer
                     ;tol)
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

function dataplot3D(xdata, ydata, zdata; plot_title,
                                         xwidth=pi/2, ywidth=pi/2)
    the_plot = plot!(xdata, ydata, zdata,
                     xlabel="square coeff", ylabel="cubic coeff",
                     zlabel="Probability",
                     seriestype=:scatter, markersize=1.5,
                     legend=:none, title=plot_title,
                     xlims=(-xwidth, xwidth), ylims=(-ywidth, ywidth),
                     zlims=(0,1),

                     )
end

"""Utility script designed to automatize plotting."""
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

"""
Copy-pasted from pss_importandrun.jl.
If you change it automatization breaks.
"""
statspath(dir, dim, keyword, uid)="$dir/random_matrices/ps$(dim)stats_$keyword$uid.csv"

function quickplot(jsonfile; kwargs...)
    model, filedata = dicttonamedtuple.(JSON.parsefile(jsonfile))
    datapath = statspath(filedata.proj_dir,
                         model.nmodes,
                         model.merit_keyword,
                         filedata.id)
    quickplot(filedata, datapath; kwargs...)
end

function quickplot3D(params, filepath; xindex=4, yindex=5, zindex=3,
                                        width=pi/2)
    data = succ_filter(readdlm(filepath, ',', skipstart=1), 2)
    succplot = dataplot3D(data[:,xindex], data[:,yindex],
                          data[:,zindex];
                          plot_title=params.plot_title,
                          xwidth=width,
                          ywidth=width
                          )
    plotpath = "$(params.plot_dir)/succ_$(params.id)_$(params.plot_title)"
    savefig(succplot, plotpath)
    println("Plot saved to \'$plotpath\'.")
end

function quickplot3D(jsonfile; kwargs...)
    model, filedata = dicttonamedtuple.(JSON.parsefile(jsonfile))
    datapath = statspath(filedata.proj_dir,
                         model.nmodes,
                         model.merit_keyword,
                         filedata.id)
    quickplot3D(filedata, datapath; kwargs...)
end

valbin(x, width, min=0) = Integer(cld(x-min, width))

function fillbins(vals, coords, shape)
    hist = Array{Union{Missing, eltype(vals)}}(missing, shape)
    for (value, indices) in zip(vals, coords)
        if ismissing(hist[indices...]) || value >= hist[indices...]
            hist[indices...] = value
        end
    end
    return hist
end

function coordlist(pointlists...; nbins, maxs, mins)
    widths = (maxs .- mins) ./ nbins
    [valbin.(i, widths, mins) for i in zip(pointlists...)]
end

"""delta is a small tolerance to prevent date form falling out of the bins
due to floating point errors"""
function binned(vals, pointlists...; nbins, delta=0.05)
    maxvec =  ntuple(x -> pi+delta, length(nbins)) #Can I not allocate an entire array?
    minvec =  ntuple(x -> -pi-delta, length(nbins))
    coords = coordlist(pointlists...; nbins=nbins, maxs=maxvec, mins=minvec)
    fillbins(vals, coords, nbins)
end

function quickheatmap(params, filepath;
                      xindex=4, yindex=5, zindex=3, devindex=2,
                      nbins=(30,30), delta=0.05, format="png"
                      )
    data = succ_filter(readdlm(filepath, ',', skipstart=1), devindex)
    zdata = binned(data[:,zindex], data[:,xindex], data[:,yindex];
                   nbins=nbins)
    succplot = heatmap(range(-pi-delta , pi+delta, length=nbins[1]),
                       range(-pi+delta, pi+delta, length=nbins[2]),
                       zdata;
                       plot_title=params.plot_title,
                       xlabel="quadratic coefficient",
                       ylabel="cubic coefficient",
                       zlabel="probability",
                       color=:greys
                       )
    plotpath="""$(params.plot_dir)/heatmap_$(params.id)_$(params.plot_title).$format"""
    display(succplot)
    savefig(succplot, plotpath)
    println("Plot saved to \'$plotpath\'.")
end

function quickheatmap(jsonfile; kwargs...)
    model, filedata = dicttonamedtuple.(JSON.parsefile(jsonfile))
    datapath = statspath(filedata.proj_dir,
                         model.nmodes,
                         model.merit_keyword,
                         filedata.id
                         )
    quickheatmap(filedata, datapath; kwargs...)
end

if !isinteractive() && (abspath(PROGRAM_FILE) == @__FILE__)
    quickplot(ARGS[1])
end
