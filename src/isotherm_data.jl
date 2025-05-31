struct AbsorbedIsothermData{TT} <: Tables.AbstractColumns
    p::Vector{TT}
    l::Vector{TT}
    T::Vector{TT}
    p_label::Symbol
    l_label::Symbol
    T_label::Symbol
end

const AdsIsoTData = AbsorbedIsothermData

#= function isotherm_data(p::P,l::L,p_label::Symbol,l_label::Symbol) where {P,L}
    return AbsorbedIsothermData(p,l,zeros(eltype(p),length(p)),p_label,l_label,:T)
end =#

function isotherm_data(p::P, l::L, p_label::Symbol, l_label::Symbol; fill_T = zeros(eltype(p))) where {P, L}
    
    return AbsorbedIsothermData(p, l, fill(fill_T, length(p)), p_label, l_label, :T) 

end

function merge_isotherm_data(tables::AdsIsoTData{TT}...) where TT
    row_tables = map(Tables.rowtable, tables)
    concatenated_rows = vcat(row_tables...)
    return Tables.columntable(concatenated_rows)
end

function isotherm_data(p::P,l::L,t::T,p_label::Symbol,l_label::Symbol,T_label::Symbol) where {P,L,T}
    @assert length(p) == length(l) == length(t)
    @assert p_label != l_label
    @assert p_label != T_label
    @assert l_label != T_label
    TT = Base.promote_eltype(p,l,t)
    n = length(p)
    pp = Vector{TT}(undef,n)
    ll = Vector{TT}(undef,n)
    tt = Vector{TT}(undef,n)

    if n > 0
        pmin = minimum(p)
        if pmin < 0
            throw(DomainError(pmin,"adsorbed input pressure must be positive"))
        end
        if !issorted(p)
            idx = sortperm(p) # Sort original p
            pp .= p[idx]
            ll .= l[idx]
            tt .= t[idx]
        else
            pp .= p
            ll .= l
            tt .= t
        end
    else # n == 0, vectors are already empty and correctly typed via undef,n. Or pp,ll,tt remain Vector{TT}(undef,0)
        # If p, l, t were Any[], TT could be Any. Forcing to Float64[] if completely empty.
        # However, if input `p` was e.g. Float64[], TT will be Float64.
        # The lines below ensure pp, ll, tt are copies, empty if p,l,t are empty.
        pp .= p
        ll .= l
        tt .= t
    end
    return AbsorbedIsothermData{TT}(pp,ll,tt,p_label,l_label,T_label)
end

isotherm_data(p::AbstractVector,l::AbstractVector) = isotherm_data(p,l,zeros(eltype(p),length(p)))
isotherm_data(p::AbstractVector,l::AbstractVector,T::AbstractVector) = isotherm_data(p,l,T,:p,:l,:T)

function isotherm_data(table,p_label::String,l_label::String,T_label::String)
    return isotherm_data(table,Symbol(p_label),Symbol(l_label),Symbol(T_label))
end

function isotherm_data(table,p_label::Symbol,l_label::Symbol,T_label::Symbol)
    if Tables.istable(table)
        p = Tables.getcolumn(table,p_label)
        l = Tables.getcolumn(table,l_label)
        T = Tables.getcolumn(table,T_label)
         return isotherm_data(p,l,T,p_label,l_label,T_label)
    else
        throw(ArgumentError("table must be a Tables.jl table"))
    end
end

function isotherm_data(table)
    return isotherm_data(table,keys(table)...)
end

#Tables.jl api, column access

Tables.istable(::Type{<:AbsorbedIsothermData}) = true
Tables.columnaccess(::Type{<:AbsorbedIsothermData}) = true
function Tables.schema(m::AbsorbedIsothermData{T}) where {T}

    Tables.Schema(Tables.columnnames(m), (T,T,T))
end
Tables.columns(m::AbsorbedIsothermData) = m

function Tables.getcolumn(table::AbsorbedIsothermData, nm::Symbol)
    p_label,l_label,T_label = getfield(table,4),getfield(table,5),getfield(table,6)
    if nm == p_label
        return getfield(table,1)
    elseif nm == l_label
        return getfield(table,2)
    elseif nm == T_label
        return getfield(table,3)
    else
        throw(KeyError(nm))
    end
end

function Tables.getcolumn(table::AbsorbedIsothermData, nm::Int)
    if nm == 1 || nm == 2 || nm == 3
        return getfield(table,nm)
    else
        throw(BoundsError(table,nm))
    end
end

function Tables.getcolumn(m::AbsorbedIsothermData, ::Type{T}, col::Int, nm::Symbol) where {T}
    return Tables.getcolumn(m, col)
end

Tables.columnnames(table::AbsorbedIsothermData) = (getfield(table,4),getfield(table,5),getfield(table,6))

Base.eltype(::Type{<:AbsorbedIsothermData{T}}) where {T} = T
Base.eltype(::AbsorbedIsothermData{T}) where {T} = T

temperature(m::AbsorbedIsothermData) = Tables.getcolumn(m,3)
loading(m::AbsorbedIsothermData) = Tables.getcolumn(m,2)
pressure(m::AbsorbedIsothermData) = Tables.getcolumn(m,1)



function split_data_by_temperature(data::AdsIsoTData{TT}) where TT
    l = loading(data)
    p = pressure(data)
    t = temperature(data)
    
    unique_temps = sort(unique(t))  # Sort unique temperatures in increasing order
    temp_data = Vector{Tuple{Vector{TT}, Vector{TT}}}(undef, length(unique_temps))
    
    for (i, temp) in enumerate(unique_temps)
        indices = findall(t .== temp)
        l_temp = l[indices]
        p_temp = p[indices]
        temp_data[i] = (l_temp, p_temp)
    end
    
    return unique_temps, temp_data
end

export AbsorbedIsothermData, isotherm_data, split_data_by_temperature, merge_isotherm_data
