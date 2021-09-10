function *!(mat::AbstractVecOrMat, bs::BeamSplitter)
    *!(bs, transpose!(mat)::AbstractVecOrMat)
    transpose!(mat)
end
