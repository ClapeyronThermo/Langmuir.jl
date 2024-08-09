Base.eltype(model::IsothermModel{T}) where T = T

function Base.eltype(::Type{M}) where M <: IsothermModel{T} where T
    return T
end

Rgas(model) = 8.31446261815324

"""
    model_length(model::IsothermModel)

Returns the amount of parameters necessary to instantiate a `model::IsothermModel`. For simple models, this is equivalent to `fieldcount(typeof(model))`, but composite models (like `MultiSite`), define their own implementation.
"""
#default.
model_length(::Type{T}) where T <: IsothermModel = _model_length(T)
model_length(model::IsothermModel) = model_length(typeof(model))

function _model_length(model::Type{T}) where T <: IsothermModel
    return fieldcount(T)
end

function isotherm_checkbounds(::Type{M},data) where M <: IsothermModel
    lb = isotherm_lower_bound(float(_eltype(data)),M)
    ub = isotherm_upper_bound(float(_eltype(data)),M)
    @inbounds for i in 1:model_length(M)
        if !(lb[i] <= data[i] <= ub[i])
            IsothermBoundsError(M,lb[i],ub[i],data[i],i)
        end
    end
    return nothing
end


function IsothermBoundsError(::Type{M},lb,ub,datai,i) where M <: IsothermModel
    description = isotherm_descriptions(M)[i]
    symbol = string(fieldname(M,i))
    if length(description) != 0
        d = lazy"($description) "
    else
        d = lazy""
    end
    throw(ArgumentError(lazy"$(nameof(M)): value for the field `$symbol` $(d)is out of the parameter bounds: ($lb <= $datai <= $ub) == false"))
end

function from_vec(m::IsothermModel,x)
    return from_vec(typeof(m),x,true)
end

from_vec(m::IsothermModel,x,check) = from_vec(typeof(m),x,check)

from_vec(::Type{M},p) where {M <: IsothermModel} = from_vec(M,p,true)

function from_vec(::Type{M},p,check) where {M <: IsothermModel}
    data = ntuple(i -> p[i], model_length(M))
    check && isotherm_checkbounds(M,data)
    return M(data...)
end

function from_vec(::Type{M},p,check) where M <: IsothermModel{T} where T
    data = ntuple(i -> T(p[i]), model_length(M))
    check && isotherm_checkbounds(M,data)
    return M(data...)
end

function to_vec!(model::IsothermModel,x)
    for i in 1:model_length(model)
        x[i] = getfield(model,i)
    end
    return x
end

function to_vec(model::IsothermModel)
    x = Vector{eltype(model)}(undef, model_length(model))
    to_vec!(model,x)
    return x
end

function to_tuple(model::IsothermModel)
    return ntuple(i -> getfield(model,i), model_length(model))
end

Base.zero(model::M) where M <: IsothermModel = Base.zero(M)

function Base.zero(model::Type{M}) where M <: IsothermModel{T} where T
    isotherm_zero(T,model)
end

#when T is not defined (zero(Langmuir))
function Base.zero(model::Type{M}) where M <: IsothermModel
    isotherm_zero(Float64,model)
end

function isotherm_zero(::Type{T},model::Type{M}) where T <: Number where M <: IsothermModel
    from_vec(M,ntuple(Returns(zero(T)),model_length(model)))
end

function Base.iszero(model::IsothermModel)
    result = true
    for i in 1:model_length(model)
        result &= iszero(getfield(model,i))
    end
    return result
end

#default: fill with ones
function x0_guess_fit(::Type{T}, data) where T <: IsothermModel
    p = pressure(data)
    _1 = one(eltype(p))
    _ones = ntuple(Returns(_1),model_length(T))
    _lb_bounds = 2 .* lower_bound(typeof(_1),T)
    v = max.(_ones,_lb_bounds)
    return from_vec(T,v)
end

function isotherm_descriptions(::Type{T}) where T <: IsothermModel
    return ntuple(Returns(""),model_length(M))
end

"""
    isotherm_lower_bound(model::IsothermModel)
    isotherm_lower_bound(T,model::IsothermModel)
    isotherm_lower_bound(T,::Type{M}) where M <:IsothermModel

Returns the lower bound for the parameters of the isotherm model `model` of type `M`. with number type `T`, as a `Ntuple{model_length(M),T}`.
The default assumes that all parameters are nonnegative.
"""
function isotherm_lower_bound(model::IsothermModel)
    return isotherm_lower_bound(typeof(model))
end

function isotherm_lower_bound(::Type{T}) where T <: IsothermModel
    return isotherm_lower_bound(_eltype(T),T)
end

function isotherm_lower_bound(::Type{T},m::IsothermModel) where T
    return isotherm_lower_bound(T,typeof(m))
end

function isotherm_lower_bound(::Type{T},::Type{M}) where T where M <: IsothermModel
    ntuple(Returns(T(-Inf)),model_length(M))
end

"""
    isotherm_upper_bound(model::IsothermModel)
    isotherm_upper_bound(T,model::IsothermModel)
    isotherm_upper_bound(T,::Type{M}) where M <:IsothermModel

Returns the upper bound for the parameters of the isotherm model `model` of type `M`. with number type `T`, as a `Ntuple{model_length(M),T}`.
The default assumes no upper bound for the parameters.
"""
function isotherm_upper_bound(model::T) where T <: IsothermModel
    return isotherm_upper_bound(typeof(model))
end

function isotherm_upper_bound(::Type{T}) where T <: IsothermModel
    return isotherm_upper_bound(_eltype(T),T)
end

function isotherm_upper_bound(::Type{T},m::IsothermModel) where T
    return isotherm_upper_bound(T,typeof(m))
end

function isotherm_upper_bound(::Type{T},::Type{M}) where T where M <: IsothermModel
    ntuple(Returns(T(Inf)),model_length(M))
end

"""
    @with_metadata(struct_expr)

macro that allows to define an isotherm model with additional metadata, about parameter bounds and descriptions of parameters:

## Usage:

```julia
AdsorbedSolutionTheory.@with_metadata struct MyIsotherm{T} <: IsothermModel{T}
    A::T,(0,1),"field A" #bounds and description provided
    B::T #nothing provided
    C::T,(1,10) #only bounds provided
    D::T,nothing,"field D" #only description provided
end

from_vec(MyIsotherm,(1,2,3,4)) #ok
from_vec(MyIsotherm,(-1,2,3,4)) #ArgumentError: MyIsotherm: value for the field `A` (field A) is out of the parameter bounds: (0.0 <= -1 <= 1.0) == false
from_vec(MyIsotherm,(1,2,-3,4)) #ArgumentError: MyIsotherm: value for the field `C` is out of the parameter bounds: (1.0 <= -3 <= 10.0) == false
```
"""
macro with_metadata(_struct_def)
    struct_def = copy(_struct_def)
    struct_fields = struct_def.args[3]
    _lb = Any[]
    _ub = Any[]

    #X{T} <: IsothermModel{T}
    struct_head_def = struct_def.args[2]
    if struct_head_def.head == :(<:)
        struct_head = struct_head_def.args[1].args[1]
    elseif struct_head_def.head == :curly
        struct_head = struct_head_def.args[1]
    else
        throw(ParseError("invalid expression for struct head: $struct_head_def"))
    end

    descriptions = String[]
    for i in 1:length(struct_fields.args)
        #it could be a LineNumberNode
        arg_i = struct_fields.args[i]

        if arg_i isa Expr
            if arg_i.head == :tuple
                tuple_expr = arg_i.args
            elseif arg_i.head == :(::)
                tuple_expr = Any[arg_i,:nothing,""] #no metadata provided
            end
            if length(tuple_expr) > 1
                bounds = tuple_expr[2] #bounds provided
            else
                bounds = :nothing #default, converted to (-Inf,Inf)
            end

            if length(tuple_expr) > 2
                description_i = tuple_expr[3] #description provided
            else
                description_i = "" #default
            end

            if bounds isa Expr && bounds.head == :tuple
                lb_i = bounds.args[1]
                ub_i = bounds.args[2]
                push!(_lb,lb_i)
                push!(_ub,ub_i)
            else
                push!(_lb,-Inf)
                push!(_ub,Inf)
            end
            push!(descriptions,description_i)
            struct_fields.args[i] = tuple_expr[1]
        end
    end
    float_lb = map(Float64,_lb)
    float_ub = map(Float64,_ub)
    lb_tuple = Expr(:tuple,float_lb...)
    ub_tuple = Expr(:tuple,float_ub...)
    descriptions_tuple = Expr(:tuple,descriptions...)
    quote
        $struct_def

        function AdsorbedSolutionTheory.isotherm_lower_bound(::Type{TT},::Type{M}) where {TT,M <: $struct_head}
            return TT.($lb_tuple)
        end

        function AdsorbedSolutionTheory.isotherm_upper_bound(::Type{TT},::Type{M}) where {TT,M <: $struct_head}
            return TT.($ub_tuple)
        end

        function AdsorbedSolutionTheory.isotherm_descriptions(::Type{M}) where {M <: $struct_head}
            return $descriptions_tuple
        end
    end |> esc
end

export isotherm_lower_bound, isotherm_upper_bound, model_length
export IsothermModel

include("freundlich.jl")
include("langmuir.jl")
include("langmuir_freundlich.jl")
include("redlich_peterson.jl")
include("sips.jl")
include("quadratic.jl")
include("bet.jl")
include("henry.jl")
include("temkin.jl")
include("unilan.jl")
include("Toth.jl")
include("multisite.jl")
include("interpolation.jl")
