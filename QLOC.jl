module QLOC

import Base: *, transpose, inv, conj, size, rand, getindex, show, Matrix

export BeamSplitter, OxfordDecomp, btoQTikz, diagonality, qtikzbody

include("Permanent.jl")
include("beamsplitter.jl")
include("reck-oxford.jl")

end

