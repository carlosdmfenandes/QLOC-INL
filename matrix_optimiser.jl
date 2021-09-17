"Optimise matrix post selection process."

using Optim

include("matrix_composer.jl")
include("csvread.jl")
include("post_selection_statistics.jl")

const HOME = homedir()
const PROJDIR = "$HOME/Documents/Quantum_Information"
const DIM = 4
const FILESDIR = "$PROJDIR/random_matrices/dim_$DIM"
const NUM = 2
"Contains a list of functions and what to name the files
this algorithm generates."
ang = pi / 4
NPHOTONS=2
NNONLIN=2
 
"Determine the name of the file to import."
csvname(dim, complex, suffix)="$FILESDIR/Unitary_2$(dim)$(complex)_example$suffix.csv"
csvname(complex, suffix)=csvname(DIM, complex, suffix)

"Determine the file path to export results."
statspath(keyword, dim=DIM, dir=PROJDIR, uid=NUM)="$dir/random_matrices/ps$(dim)stats_$keyword$uid.csv"

"Find the directory matrix to export."
matdir(keyword, dim=DIM, dir=PROJDIR, uid=NUM)="$dir/random_matrices/ps$(dim)_$keyword$uid"

"Wrapper around optimize to convert matrices to parametrized form and back"
function ps_optimizer(init_matrix, merit_function, function_args=(), optimizeargs...)
    inits = argument_form(init_matrix)
    result = optimize(inits, optimizeargs...) do x
         merit_function(nsamplitudes(matrix_form(x), NNONLIN),
                        function_args...)
    end
    return result
end
function ps_optimizer(init_matrix, keyword::AbstractString, args...) 
    ps_optimizer(init_matrix, MERIT_DICTS[keyword], args...) 
end

"""
Extract the physically relevant results from an Optim.Optimization object 
returned by ps_optimizer. 
Return the minimal matrix and a vector containing the deviation form uniformity, 
the non-linear angle and the intererometer's success probability.
"""
function ps_results(optimization)
    min = Optim.minimizer(optimization)
    minmat = matrix_form(min, DIM)
    amplitudevec = nsamplitudes(minmat, NNONLIN) 
    mean_amplitude = mean(abs.(amplitudevec))
    deviation = std(amplitudevec, mean=mean_amplitude)
    rel_dev = deviation / mean_amplitude
    nonlinarg = anglestocoefs(angle.(amplitudevec))[NNONLIN]
    succ_prob = abs2(mean_amplitude)
    minmat, [rel_dev; nonlinarg; succ_prob]
end

"""
Solve an optimization problem and write the results to a file. 
Takes as input the csv files containing the real and imaginary parts of the input matrix and 
a function to_optimize. Alternatively this function can be specified by an
apropriate keyword as defined in the MERIT_DICTS dictionary in the \"post_selection_statistics.jl\" file.
It can also pass additional arguments to the function to_optimize.
The files to which the results are to be written must be specified by the stats_file and matrix_file 
keywords.
"""   
function writetofile(real_path, imag_path, to_optimize,  args...; stats_file=devnull, matrix_file=devnull, prefix=nothing)
    init_matrix = matread(real_path, imag_path)
    optimres = ps_optimizer(init_matrix, to_optimize, args...)
    minmat, data = ps_results(optimres)
    if isnothing(prefix)
        writedlm(stats_file, data, ',')
    else 
        writedlm(stats_file, [prefix; data], ',')
        println(matrix_file, prefix)
    end
    writedlm(matrix_file, minmat, ',')
end

function loop(keyword="stddev", some_range=1001:100000)
    statsfile = statspath(keyword)
    filestream = open(statsfile, "w")
    try
        for x in some_range
            real_matrix_readpath = csvname("R",x)
            imag_matrix_readpath = csvname("I",x)
	    mdir = matdir(keyword)
            matrix_filepath = "$(mkpath(mdir))/$x.csv"
            matrix_filestream = open(matrix_filepath, "w+")
            try
                writetofile(real_matrix_readpath,
                            imag_matrix_readpath,
                            keyword,
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

if !isinteractive()
    loop()
end

