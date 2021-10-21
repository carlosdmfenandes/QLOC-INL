"Optimise matrix post selection process."

using Optim

include("matrix_composer.jl")
include("csvread.jl")
include("post_selection_statistics.jl")

"Wrapper around optimize to convert matrices to parametrized form and back"
function ps_optimizer(init_matrix, merit_function, optimargs...; 
                      nphotons, function_args=NamedTuple(), optimkwargs...)
    inits = argument_form(init_matrix)
    result = optimize(inits, optimargs...; optimkwargs...) do x
        merit_function(nsamplitudes(matrix_form(x), nphotons);
                       function_args...)
        end
    return result
end
function ps_optimizer(init_matrix, keyword::AbstractString, args...; kwargs...) 
    ps_optimizer(init_matrix, MERIT_DICTS[keyword], args...; kwargs...) 
end

"""
Extract the physically relevant results from an Optim.Optimization object 
returned by ps_optimizer. 
Return the minimal matrix and a vector containing the deviation form uniformity, 
the non-linear angle and the interferometer's success probability.
"""
function ps_results(optimization; nphotons)
    min = Optim.minimizer(optimization)
    minmat = matrix_form(min) #Use explicit dimension if it breaks.
    amplitudevec = nsamplitudes(minmat, nphotons) 
    mean_amplitude = mean(abs.(amplitudevec))
    deviation = std(abs.(amplitudevec), corrected=false, mean=mean_amplitude)
    rel_dev = deviation / mean_amplitude
    nonlinarg = anglestocoeffs(angle.(amplitudevec))[3:(nphotons+1)]
    succ_prob = abs2(mean_amplitude)
    minmat, [rel_dev; succ_prob; nonlinarg]
end
