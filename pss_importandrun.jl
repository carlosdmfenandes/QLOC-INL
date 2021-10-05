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
a function to_optimize. Alternatively this function can be specified by an
apropriate keyword as defined in the MERIT_DICTS dictionary in the \"post_selection_statistics.jl\" file.
It can also pass additional arguments to the function to_optimize.
The files to which the results are to be written must be specified by the stats_file and matrix_file 
keywords.
"""   
function writetofile(real_path, imag_path, to_optimize,  args...; 
                     stats_file=devnull, matrix_file=devnull,
                     prefix=nothing, nphotons=2, fargs=())
    init_matrix = matread(real_path, imag_path)
    optimres = ps_optimizer(init_matrix, to_optimize, args...; 
                            function_args=fargs)
    minmat, data = ps_results(optimres)
    if isnothing(prefix)
        writedlm(stats_file, reshape(data,1,:), ',')
    else 
        writedlm(stats_file, reshape(Any[prefix;data],1,:), ',') #data would be better as tuple. Rethink?
        println(matrix_file, prefix)
    end
    writedlm(matrix_file, minmat, ',')
end

function calc(opt_params, file_params, some_range)
    statsfile = statspath(file_params.proj_dir,
                opt_params.nmodes,
                opt_params.merit_keyword,
                file_params.id)
    filestream = open(statsfile, "w+")
    csvdir= "$(file_params.proj_dir)/random_matrices/dim_$(opt_params.nmodes)" 
    header = ["index" "relative deviation" "success probability" "nonlinearity"]
    writedlm(filestream, header, ',')
    try
        for x in some_range
            real_matrix_readpath = csvname(file_params.import_dir, opt_params.nmodes, "R",x)
            imag_matrix_readpath = csvname(file_params.import_dir, opt_params.nmodes, "I",x)
            mwdir = matwritedir(file_params.proj_dir, opt_params.nmodes, 
                          opt_params.merit_keyword, file_params.id)
            matrix_filepath = "$(mkpath(mwdir))/$x.csv"
            matrix_filestream = open(matrix_filepath, "w+")
            try
                writetofile(real_matrix_readpath,
                            imag_matrix_readpath,
                            opt_params.merit_keyword;
                            fargs = get(opt_params, :args, ()),
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

#=This function is type unstable but it doesn't seem to have a 
performance penalty compared to bare calc.
=#
function importandrun(paramsfile, irange)
    paramdata = dicttonamedtuple.(JSON.parsefile(paramsfile))
    (optpar, filespar) = (paramdata[1], paramdata[2]) 
    calc(optpar, filespar, irange)
end

#=

Originally intended to correct type instabilities in importandrun, 
but there seem not to be a performace penalty.

include("parameter_type.jl")

function importandrunstable(paramsfile, irange, args...)
    paramdata = JSON.parsefile(paramsfile)
    optpar = PssModel(; symzip(paramdata[1])... )
    filespar = PssFiles(; symzip(paramdata[2])... ) 
    calc(optpar, filespar, irange)
end

=#
