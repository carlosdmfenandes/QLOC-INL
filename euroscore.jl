const POINTS=[12, 10, 8, 7, 6, 5, 4, 3, 2, 1]

function addscore!(score, vote)
    for (i,x) in enumerate(vote)
        score[x] = get(score, x, 0) + POINTS[i]
    end
    score
end

function scorecount(votes)
    score = Dict{eltype(votes[1]), Int64}()
    for i in votes
        addscore!(score, i)
    end
    score
end

#Test code
#
vote1 = ['q', 'w', 'e', 'r', 't', 'y', 'u' , 'i' ,'o', 'p']
vote2 = ['z', 'a', 'q', 'w', 's', 'x', 'c', 'd', 'e', 'r']
vote3 = ['a', 'q', 'w', 's', 'd', 'e', 'r', 'f', 'g', 't']
tvote = [vote1, vote2, vote3]
#
