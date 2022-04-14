module QLOC

import Base: *, transpose, inv, conj, size, rand

export BeamSplitter

include("permanent.jl")
include("pairwise.jl")
include("beamsplitter.jl")
end

