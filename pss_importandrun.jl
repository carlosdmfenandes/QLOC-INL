"Imports data and calculates the model." 

import JSON

include("matrix_optimiser.jl")

"Determine the name of the file to import."
csvname(dir, dim, complex, suffix)="$dir/Unitary_2$(dim)$(complex)_example$suffix.csv"

"Determine the file path to export results."
statspath(dir, dim, keyword, uid)="$dir/random_matrices/ps$(dim)stats_$keyword$uid.csv"

"Find the directory matrix to export."
matwritedir(dir, dim, keyword, uid)="$dir/random_matrices/ps$(dim)_$keyword$uid"

"Covert a dictionary to a Named Tuple"
dicttonamedtuple(dic::AbstractDict) = (; zip(Symbol.(keys(dic)),values(dic))...)

symzip(dic)=zip(Symbol.(keys(dic)),values(dic))

"""
Solve an optimization problem and write the results to a file. 
Takes as input the csv files containing the real and imaginary parts of the input matrix and 
a function to_optimize. 
Alternatively this function can be specified by an
apropriate keyword as defined in the MERIT_DICTS dictionary in the \"post_selection_statistics.jl\" file.
It can also pass additional arguments to the function to_optimize.
The files to which the results are to be written must be specified by the stats_file and matrix_file 
keywords.
"""   
function writetofile(real_path, imag_path, to_optimize,  args=(); 
                     stats_file=devnull, matrix_file=devnull,
                     prefix=nothing, nphotons, fargs=(;) )
    init_matrix = matread(real_path, imag_path)
    optimres = ps_optimizer(init_matrix, to_optimize, args...; 
                            nphotons=nphotons, function_args=fargs)
    minmat, data = ps_results(optimres; nphotons=nphotons)
    if isnothing(prefix)
        writedlm(stats_file, reshape(data,1,:), ',')
    else 
        writedlm(stats_file, reshape(Any[prefix;data],1,:), ',') #data would be better as tuple. Rethink?
        println(matrix_file, prefix)
    end
    writedlm(matrix_file, minmat, ',')
end

function solvewtf(init_matrix, to_optimize,  args=(); 
                     stats_file=devnull, matrix_file=devnull,
                     prefix=nothing, nphotons, fargs=(;) )
    optimres = ps_optimizer(init_matrix, to_optimize, args...; 
                            nphotons=nphotons, function_args=fargs)
    minmat, data = ps_results(optimres; nphotons=nphotons)
    if isnothing(prefix)
        writedlm(stats_file, reshape(data,1,:), ',')
    else 
        writedlm(stats_file, reshape(Any[prefix;data],1,:), ',') #data would be better as tuple. Rethink?
        println(matrix_file, prefix)
    end
    writedlm(matrix_file, minmat, ',')
end

"""
Takes the optimization parameters and reads the 
arguments to be passed to the merit function.
"""
function arg_manager(opt_params)
    keyword=opt_params.merit_keyword
    if keyword == "setpar"  
        return (opt_params.setangs, opt_params.weights)
    else
        return ()
    end
end

"""
Takes named tuples or structs of parameters and solves the 
optimization sampling problem. some_range specifies the 
range of matrices stored as files to use.
"""
#IT'S A FUCKING MESS
function calc(opt_params, file_params, some_range; wmode="w+", fargs=(;) )
    statsfile = statspath(
                file_params.proj_dir,
                opt_params.nmodes,
                opt_params.merit_keyword,
                file_params.id
                )
    filestream = open(statsfile, wmode)
    csvdir= "$(file_params.proj_dir)/random_matrices/dim_$(opt_params.nmodes)" 
    header = ["index" "relative deviation" "success probability" "nonlinearity"]
    writedlm(filestream, header, ',')
    try
        for x in some_range
            real_matrix_readpath = csvname(file_params.import_dir,
                                           opt_params.nmodes, "R",x)
            imag_matrix_readpath = csvname(file_params.import_dir,
                                           opt_params.nmodes, "I",x)
            mwdir = matwritedir(file_params.proj_dir, opt_params.nmodes, 
                               opt_params.merit_keyword, file_params.id)
            matrix_filepath = "$(mkpath(mwdir))/$x.csv"
            matrix_filestream = open(matrix_filepath, "w+")
            init_matrix=matread(real_path, imag_path)
            try
                solvewtf(init_matrix,
                        opt_params.merit_keyword;
                        fargs = fargs,
                        nphotons = opt_params.nonlinphots,
                        stats_file=filestream,
                        matrix_file=matrix_filepath,
                        prefix=x
                        )
            finally
                close(matrix_filestream)
            end
        end
    finally
        close(filestream)
    end
    print("done")
    return 0
end

function multicalc(opt_params, file_params, matrixfiles; 
                   fargs=(;))
    statsfile = statspath(
                file_params.proj_dir,
                opt_params.nmodes,
                opt_params.merit_keyword,
                file_params.id
                )
    if isfile(statsfile)
        filestream = open(statsfile, "a") #This function is to be looped.
    elseif ispath(statsfile)
        error("Cannot write to $statsfile. Path exits, but is not a file.")
    else
        filestream = open(statsfile, "a+") #This function is to be looped.
        header = ["relative deviation" "success probability" "quadratic" "cubic"]
        writedlm(filestream, header, ',')
    end
    csvdir= "$(file_params.proj_dir)/random_matrices/dim_$(opt_params.nmodes)" 
    try
        for init_file in matrixfiles
            init=readdlm(init_file, ',', ComplexF64)
            solvewtf(init,
                    opt_params.merit_keyword;
                    fargs = fargs,
                    nphotons = opt_params.nonlinphots,
                    stats_file=filestream,
                    )
        end
    finally
    close(filestream)
    end
    println("Done writing to $statsfile.")
    return 0
end

#=
This function is type unstable but it doesn't seem to have a 
performance penalty compared to bare calc.
=#
function importandrun(paramsfile, irange)
    jsondicts = JSON.parsefile(paramsfile)
    if length(jsondicts) == 3
        fargs = convert(Dict{String, Vector{Int}} ,jsondicts[3])
        fargs = dicttonamedtuple(fargs)
    else
        fargs = (;)
    end
    paramdata = dicttonamedtuple.(jsondicts)
    (optpar, filespar) = (paramdata[1], paramdata[2])
    calc(optpar, filespar, irange; fargs=fargs)
end

function importandrun(paramsfile, irange, fargs=(;))
    jsondicts = JSON.parsefile(paramsfile)
    paramdata = dicttonamedtuple.(jsondicts)
    (optpar, filespar) = (paramdata[1], paramdata[2])
    calc(optpar, filespar, irange; fargs=fargs)
end

"""Generate an iterator spanning all files (excluding directories)
in dir and its subdirectories."""
iterfiles(dir) =   (joinpath(i,f) 
                        for (i,j,files)=walkdir(dir) 
                            for f=files
                    )
                               
function importtohist(paramsfile, matdir; 
                      nbins, trials)
    jsondicts = JSON.parsefile(paramsfile)
    paramdata = dicttonamedtuple.(jsondicts)
    (optpar, filespar) = paramdata[1:2]
    length(paramdata) >=3 ? kpars=paramdata[3] : kpars=(;)
    rangevals = range(-pi, pi, length=nbins)
    matrices =  Iterators.take(iterfiles(matdir), trials)
    for (a2, a3) in  Iterators.product(rangevals, rangevals)
        multicalc(optpar, filespar, matrices; 
                  fargs=(;weights=[0., 0, 1, 1], 
                         phases=[0, 0, a2, a3], 
                         kpars...)
                  )
    end
end


if !isinteractive() && (abspath(PROGRAM_FILE) == @__FILE__)
    importandrun(ARGS[1], ARGS[2])
end

#=
Originally intended to correct type instabilities in importandrun, 
but there seems not to be a difference in performace .

include("parameter_type.jl")

function importandrunstable(paramsfile, irange, args...)
    paramdata = JSON.parsefile(paramsfile)
    optpar = PssModel(; symzip(paramdata[1])... )
    filespar = PssFiles(; symzip(paramdata[2])... ) 
    calc(optpar, filespar, irange)
end

=#
