include("QLOC.jl")

using .QLOC

const FAIL_CODE = 1

function event_prob(
    circuit::AbstractMatrix;
    istate=nothing, ostate=nothing, smatrix=nothing
    )
    if istate |> isnothing
        if ostate |> isnothing
            result = abs2(permanent(circuit))
        else
            result = abs2(multipermanent(circuit, istate))/prod(factorial, ostate)
        end
    else
        if ostate |> nothing
            exit(FAIL_CODE)
        else
            result = abs2(multipermanent(circuit, istate, ostate))/prod(factorial, ostate)/prod(factorial, istate)
        end
    end
end
