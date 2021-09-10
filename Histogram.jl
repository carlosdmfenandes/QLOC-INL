using DelimitedFiles
using CSV
using LinearAlgebra
using StatsPlots

sourcedir = "/home/carlosfernandes/Documents/Quantum_Information/random_matrices"
path1 = joinpath(sourcedir,"Permanents2_dim3_subdim3.csv")
path2 = joinpath(sourcedir,"Permanents2_dim4_subdim3.csv")
path3 = joinpath(sourcedir,"Permanents2_dim5_subdim3.csv")
path4 = joinpath(sourcedir,"Permanents2_dim6_subdim3.csv")
path1_ = joinpath(sourcedir,"Permanents2_dim7_subdim3.csv")
path2_ = joinpath(sourcedir,"Permanents2_dim8_subdim3.csv")
path3_ = joinpath(sourcedir,"Permanents2_dim9_subdim3.csv")
path4_ = joinpath(sourcedir,"Permanents2_dim10_subdim3.csv")
data1 = Float64.(readdlm(path1,',',  skipstart=1)[:,3])
data2 = Float64.(readdlm(path2,',',  skipstart=1)[:,3])
data3 = Float64.(readdlm(path3,',',  skipstart=1)[:,3])
data4 = Float64.(readdlm(path4,',',  skipstart=1)[:,3])
data5 = Float64.(readdlm(path1_,',',  skipstart=1)[:,3])
data6 = Float64.(readdlm(path2_,',',  skipstart=1)[:,3])
data7 = Float64.(readdlm(path3_,',',  skipstart=1)[:,3])
data8 = Float64.(readdlm(path4_,',',  skipstart=1)[:,3])
plot_data = [data1 data2 data3 data4 data5 data6 data7 data8]
data_labels = ["n=3" "n=4" "n=5" "n=6" "n=7" "n=8" "n=9" "n=10"]

fig=plot(plot_data, label=data_labels, fillalpha =1,
        normed=:probability, title="3 Photons",
         xguide="Gap", yguide="Fraction", size=(900,600),
        yscale=:log10, ylims=(0.0001,0.1), seriestype=:stephist)

export_name =  "4ph_Kr_hist.png"
path2 = joinpath(sourcedir,"ratiosTab_ph4_dim4.csv")
path3 = joinpath(sourcedir,"ratiosTab_ph4_dim5.csv")
path4 = joinpath(sourcedir,"ratiosTab_ph4_dim6.csv")
path2_ = joinpath(sourcedir,"quantumK_Tab_ph4_dim4.csv")
path3_ = joinpath(sourcedir,"quantumK_Tab_ph4_dim5.csv")
path4_ = joinpath(sourcedir,"quantumK_Tab_ph4_dim6.csv")
data2 = readdlm(path2_,',')/4#./readdlm(path2,',')
data3 = readdlm(path3_,',')/4#./readdlm(path3,',')
data4 = readdlm(path4_,',')/4#./readdlm(path4,',')
plot_data = [data2 data3 data4]
data_labels = ["n=4" "n=5" "n=6"]

plot(plot_data, label=data_labels, fillalpha=1, linewidth=1,
     seriestype=:density, bar_position="stack")

savefig(fig, "$sourcedir/3phgapplot.png")
