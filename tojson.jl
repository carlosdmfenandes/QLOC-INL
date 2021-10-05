import JSON

if isinteractive()
pfile = readline()
else
pfile = ARGS[1]
end

include(pfile)

text = JSON.json([OPTPARAMS, FILEPARAMS], 1)
path = "$(splitext(pfile)[1]).json"
io = open(path, "w+")

if abspath(PROGRAM_FILE) == @__FILE__
    try
        write(io, text)
    finally
        close(io)
    end
end
