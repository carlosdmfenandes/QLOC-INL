using StatsPlots
using SpecialFunctions

include("Permanent.jl")

cosmat(th) = [cos(th) sin(th); -sin(th) cos(th)]

gfac(x) = gamma(x+1)

function eventprob(m, xvec, yvec)
    num = abs2(rmultipermanent(m, xvec, yvec))
    den1 = prod(gfac.(xvec))
    den2 = prod(gfac.(yvec))
    return (num/den1)/den2
end

"Overflow resistant version of eventprob"
function eventprob_(m, xvec, yvec; B=true)
    tot = sum(xvec)
    num = abs2(rmultipermanent(m*ℯ/tot, xvec, yvec; B=B))
    imed = 2*log(tot/ℯ)*tot-sum(logfactorial.(xvec))-sum(logfactorial.(yvec))
    fac = exp(imed)
    return num*fac
end

a = cosmat(rand(BigFloat))
a = cosmat(BigFloat(pi)/4)
dataarray = [eventprob_(a, [100-i; i], [50; 50]) for i=0:100]
plot(0:100, dataarray, seriestype=:sticks, legend=:none)
dataarray_ = [eventprob_(a, [100-i; i], [50; 50], B=true) for i=0:100]
plot(0:100, dataarray_, seriestype=:sticks, legend=:none, lc=:red)

angs = BigFloat[0, 1//11, 1//4, 1//2, 3//4, 10//11]./2
for th in angs
    mat = cosmat(pi*th)
    dataarray_ = [eventprob_(mat, [100-i; i], [100; 0], B=true) for i=0:100]
    pl = plot(0:100, dataarray_, seriestype=:sticks, title=string(th), lc=:green)
    filename = "sshtransport_"*string(th)[1:min(end,4)]*"_.png"
    savefig(pl, filename)
end
