const LINESEP = "\\\\"

gaterep(size, label) = "\\gate[$(size)]{$(label)}"

function joinlines(mstr)
    [join(push!(strar, LINESEP), " & ", raw" & \qw ") for strar in mstr]
end

"
Convert an array of beam splitters to code for the QTikz LaTeX package
that draws a diagram of the array.
"
function toQTikz(body::Vector{String})
    preamble = [
        raw"\documentclass[tikz]{standalone}",
        raw"\usepackage{tikz}",
        raw"\usetikzlibrary{quantikz}",
    ]
    beginwrapper = [
        raw"\begin{document}",
        raw"\begin{quantikz}[transparent]",
    ]
    endwrapper = [
        raw"\end{quantikz}",
        raw"\end{document}"
    ]
    join(append!(preamble, beginwrapper, body, endwrapper), '\n')
end

toQTikz(x) = toQTikz(qtikzbody(x))

"
Return the body of a qtikz block representing an array `bsvec` of
BeamSplitter objects.
"
function qtikzbody(bsvec, header=[])::Vector{String}
    nlines = maximum(bs->max(bs.n, bs.m), bsvec)
    mstr = [String[] for i=1:nlines]
    for (l, h) in zip(mstr, header)
        append!(l,h)
    end
    for bs in bsvec
        for (index, array) in pairs(mstr)
            push!(array, qtikzcell(bs, index))
        end
    end
    joinlines(mstr)
end

#=
function qtikzcell(bs::BeamSplitter, row_num::Int)
    low, hi = minmax(bs.n, bs.m)
    anglestring = @sprintf("%4.2f", bs.rangle)
    if row_num == low
        str = gaterep(hi-low+1, anglestring)
    elseif low < row_num < hi
        str = raw"\strikethrough"
    else
        str = raw"\qw"
    end
    str
end
=#

function qtikzph(bs::BeamSplitter, row_num::Int)
    if row_num == bs.m
        phstring = @sprintf("%4.2f", bs.phase)
        phrep = gaterep(1, phstring)
    else
        phrep = raw"\qw"
    end
end

"""
Represent the Beamsplitter. A bit cryptic for the case where
bs.m > bs.n. The only difference is the mode in which the phase
shifter acts.
"""
function qtikzsplit(bs::BeamSplitter, row_num::Int)
    low, hi = minmax(bs.m, bs.n)
    if row_num == low
        reg_angle = bs.m == low ? bs.rangle : -bs.rangle
        anglestring = @sprintf("%4.2f", reg_angle)
        anglerep = gaterep(hi-low+1, anglestring)
    elseif low < row_num < hi
        anglerep = raw"\strikethrough"
    else
        anglerep = raw"\qw"
    end
end

"
Find the elements of a column representing
the action of a single beam splitter.
"
function qtikzcell(bs::BeamSplitter, row_num::Int)
    phrep = qtikzph(bs::BeamSplitter, row_num::Int)
    anglerep = qtikzsplit(bs::BeamSplitter, row_num::Int)
    bs.transposed ? "$phrep & $anglerep" : "$anglerep & $phrep"
end

qtikzphaseline(phases) = qtikzphaseline(phases, length(phases))

function qtikzphaseline(phases, nlines::Int)
    mstr = [String[] for i=1:nlines]
    for (ph, strvec) in zip(phases, mstr)
        phstring = @sprintf("%4.2f", ph)
        phgate = gaterep(1, phstring)
        push!(strvec, phgate)
    end
    return mstr
end

"
Return the body of a qtikz block representing an array `bsvec` of
`BeamSplitter` objects.
"
function qtikzbody(od::OxfordDecomp)
    len = length(od.diag)
    mstr = qtikzphaseline(angle.(od.diag), len)
    for bs in Iterators.reverse(od.bsvec)
        for (i, strvec) in pairs(mstr)
            if bs.transposed
                pushfirst!(strvec, qtikzcell(bs, i))
            else
                push!(strvec, qtikzcell(bs, i))
            end
        end
    end
    joinlines(mstr)
end

