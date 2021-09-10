newton_iter(f, df, x0) = x0 - f(x0)/df(x0)
function newton_method(f, df, x0, tol, max_iters=100)
    xn = x0
    iters = 0
    next = max_iters > 0
    while next
        xm, xn = xn, newton_iter(f, df, xn)
        iters += 1
        next = (iters < max_iters) && (abs(xn-xm)/(1+abs(xn)) >= tol)
    end
    return xn, iters
end
