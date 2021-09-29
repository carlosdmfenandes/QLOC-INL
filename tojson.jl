import JSON

pfile ="pss_parameter_files/dim3-ph2-j.jl"

include(pfile)

text = JSON.json([OPTPARAMS, FILEPARAMS])
path = "$(splitext(pfile)[1]).json"
io = open(path, "w+")
try
    write(io, text)
finally
    close(io)
end
