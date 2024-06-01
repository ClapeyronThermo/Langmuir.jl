function halley_fixpoint(f,x)
    f,df,f2f = f∂f∂2f(f,x)
    return x - 2*f*df/(2*df*df - f*d2f)
end

include("nlsolve.jl") #api for nlsolvers.jl
include("iast.jl") #fastIAS
include("fit.jl") #fitting isotherms to data
include("reverse_iast.jl") #reverse IAST (#TOD)

