module QLOC

import Base: *, transpose, inv, conj, size, rand, getindex, show, Matrix

export BeamSplitter, OxfordDecomp, toQTikz, diagonality

include("Permanent.jl")
include("beamsplitter.jl")
include("reck-oxford.jl")
include("QLOCtoTikz.jl")

end

