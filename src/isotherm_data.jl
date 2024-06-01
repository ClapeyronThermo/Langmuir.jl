struct AbsorbedIsothermData{T} <: Tables.AbstractColumns
    p::Vector{T}
    l::Vector{T}
    x_label::Symbol
    y_lable::Symbol
end

const AdsIsoTData = AbsorbedIsothermData

function isotherm_data(x::X,y::Y,x_label::Symbol,y_label::Symbol) where {X,Y}
    @assert length(x) == length(y)
    @assert x_label != y_label
    pmin = minimum(x)
    if pmin < 0
        throw(DomainError(pmin,"adsorbed input pressure must be positive"))
    end
    T = Base.promote_eltype(x,y)
    n = length(x)
    xx = Vector{T}(undef,n)
    yy = Vector{T}(undef,n)
    
    if !issorted(x)
        idx = sortperm(xx)
        yv,xv = @view(y,idx),@view(x,idx)
        yy .= yv
        xx .= xv
    else
        xx .= x
        yy .= y
    end

    return AbsorbedIsothermData{T}(x,y,x_label,y_label)
end

isotherm_data(x::AbstractVector,y::AbstractVector) = isotherm_data(x,y,:p,:l)

function isotherm_data(table,x_label::String,y_label::String)
    return isotherm_data(x,y,Symbol(x_label),Symbol(y_label))
end

function isotherm_data(table,x_label::Symbol,y_label::Symbol)
    if Tables.istable(table)
    x = Tables.get_column(table,x_label)
    y = Tables.get_column(table,y_label)
    return isotherm_data(x,y,x_label,y_label)
    else
        throw(ArgumentError("table must be a Tables.jl table"))
    end
end

#Tables.jl api, column access

Tables.istable(::Type{<:AbsorbedIsothermData}) = true
Tables.columnaccess(::Type{<:AbsorbedIsothermData}) = true
Tables.schema(m::AbsorbedIsothermData{T}) where {T} = Tables.Schema([m.x_label,m.y_label], [T,T])
Tables.columns(m::AbsorbedIsothermData) = m

function Tables.getcolumn(table::AbsorbedIsothermData, nm::Symbol)
    x_label,y_label = table.x_label,table.y_label
    if nm == x_label
        return table.x
    elseif nm == y_label
        return table.y
    else
        throw(KeyError(nm))
    end
end

function Tables.getcolumn(table::AbsorbedIsothermData, nm::Int)
    if nm == 1 || nm == 2
        return getfield(table,nm)
    else
        throw(BoundsError(table,nm))
    end
end

Tables.getcolumn(m::AbsorbedIsothermData, ::Type{T}, col::Int, nm::Symbol) where {T} = Tables.getcolumn(m,nm)


Tables.columnnames(table) = (table.x_label,table.y_label)

export AbsorbedIsothermData, isotherm_data