using DelimitedFiles
using LinearAlgebra
using StatsPlots

classical(n) = prod(i//n for i=1:n-1)

datax = [2 4 8 16]
datay = [1/2 -5/32 -837/131072 -77782159637/562949953421312]
datac = [classical(i) for i in datax]
datar = 1 ./(1 .-(datay./datac))
p = plot(datax, datar, seriestype=:scatter, xlabel="dims", leg=:none)
savefig("hadamard_ratio_gap.svg")
p
