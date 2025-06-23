function _format_param_value(value::AbstractFloat, sigdigits::Int)
    if iszero(value)
        return "0.0" # Consistent display for zero
    else
        # %.Xg formats to X significant digits, choosing %f or %e automatically
        return Printf.@sprintf("%.*g", sigdigits, value)
    end
end

# Handle other number types (like Integers) simply
function _format_param_value(value::Number, sigdigits::Int) # sigdigits is ignored here
    return string(value)
end

# Handle non-numeric types using repr (e.g., for Strings, Symbols if they occur)
function _format_param_value(value, sigdigits::Int) # sigdigits is ignored here
    return repr(value)
end


function Base.show(io::IO, a::MIME"text/plain", model::MultiSite)

    isotherms = model.isotherms
    #print(io, "$model")   

    for (idx, isotherm) in enumerate(isotherms)
        print(io, "Site $idx: ")
        Base.show(io, a, isotherm)
        ifelse(idx < length(isotherms), print(io, "\n\n"), continue)
    end
end

function Base.show(io::IO, ::MIME"text/plain", model::IsothermModel{T}) where {T}
    # 1. Print the type name and its parameter type
    # Using T from the signature is cleaner than typeof(model).parameters[1]
    print(io, "$model")

    typeofmodel = typeof(model)

    description = isotherm_descriptions(typeofmodel)

    fields = fieldnames(typeofmodel)

    for (idx, fname) in enumerate(fields)

        value = getfield(model, fname)
        desc = description[idx]
        
        value_str = _format_param_value(value, 6)

        print(io, "\n ", fname, " (", desc, "): ", value_str)

    end
     
end

function Base.show(
    io::IO,
    ::MIME"text/plain",
    data::AbsorbedIsothermData;
    allrows=!get(io, :limit, false),
    allcols=!get(io, :limit, false),
)
    nrow = length(getfield(data, 1))
    
    # Show summary line
    print(io, "AbsorbedIsothermData with $(nrow) data points")
    
    nrow == 0 && return nothing
    
    println(io)
    
    if allcols && allrows
        crop = :none
    elseif allcols
        crop = :vertical
    elseif allrows
        crop = :horizontal
    else
        crop = :both
    end
    
    # Use data directly as a Tables.jl instance
    return pretty_table(
        io,
        data;  # Use data directly since it implements Tables.jl interface
        newline_at_end=false,
        reserved_display_lines=2,
        header_alignment=:l,
        crop=crop,
        vcrop_mode=:middle,
        formatters=(v,i,j) -> 
            j == 1 ? @sprintf("%.3f", v) :  # Pressure
            j == 2 ? @sprintf("%.3f", v) :  # Loading
            j == 3 ? @sprintf("%.2f", v) :  # Temperature
            j == 4 ? @sprintf("%.3f", v) :  # Uncertainty
            v
    )
end