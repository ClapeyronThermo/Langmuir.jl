"""
    plot(model::IsothermModel, T::Number, p_range::Tuple{<:Number, <:Number}; npoints = 100)

    Plots the isotherm loading as a function of pressure for a given model and temperature.
    The function generates a plot with the specified number of points between the given pressure range.

"""

@recipe function f(model::IsothermModel, T::Number, p_range::Tuple{<:Number, <:Number}; npoints = 100)
    # --- Scientific Plot Defaults ---
    # Plot Size and Resolution
    width = 500 # width in pixels
    size --> (width, width/1.618) # width, height in pixels
    dpi --> 300

    # Font Sizes
    titlefontsize --> 12
    guidefontsize --> 10 # For x/y labels
    tickfontsize --> 8
    legendfontsize --> 8

    # Frame and Grid
    framestyle --> :box # Adds a box around the plot
    grid --> true
    gridalpha --> 0.2 # Make grid lines less prominent
    minorgrid --> true
    minorgridalpha --> 0.1

    # Ticks
    tick_direction --> :in # Ticks point inwards
    minorticks --> 5 # Number of minor ticks between major ticks

    # Line Width for the plotted data
    linewidth --> 1.5 # Default line width for series from this recipe

    # Legend
    legend --> :best # Or :topright, :bottomright, etc.

    # Margins (optional, adjust as needed)
    # left_margin --> 5Plots.mm
    # bottom_margin --> 5Plots.mm
    # --- End Scientific Plot Defaults ---

    # Plot-specific attributes (can override defaults)
    # title --> "Isotherm Plot" # Generic title, or set by a calling recipe
    xlabel --> "Pressure [Pa]"
    ylabel --> "Loading [mol/kg]"
    
    # Calculate pressure range and loadings
    p_min, p_max = p_range
    ps = range(p_min, stop=p_max, length = npoints)
    loadings = map(p_val -> loading(model, p_val, T), ps)
    
    # Label for the current series: Model Name and Temperature
    model_name = typeof(model).name.name
    label --> Printf.@sprintf("%s, T = %.2f K", model_name, T) 

    @series begin
        
    # Return x and y data for the series
    seriestype := :line
    ps, loadings
    end 

    return nothing

    end


@recipe function f(data::AdsIsoTData, T::Number) # T is now a mandatory Number
    # --- Scientific Plot Defaults ---
    width = 500
    size --> (width, Int(round(width/ 1.618))) # Golden ratio
    dpi --> 300
    titlefontsize --> 12
    guidefontsize --> 10
    tickfontsize --> 8
    legendfontsize --> 8
    framestyle --> :box
    grid --> true
    gridalpha --> 0.2
    minorgrid --> true
    minorgridalpha --> 0.1
    tick_direction --> :in
    minorticks --> 5
    # --- End Scientific Plot Defaults ---

    # Plot-specific attributes
    xlabel --> "Pressure [Pa]" # Use stored label
    ylabel --> "Loading [kg/mol]"   # Use stored label
    legend --> :best

    # Extract all data and split by temperature
    # Assumes split_data_by_temperature is accessible
    unique_temps_all, temp_specific_data_all = split_data_by_temperature(data)

    # Find the index of the target temperature
    # Allow for small floating point inaccuracies
    idx = findfirst(t_val -> isapprox(t_val, T, atol=1e-2), unique_temps_all) 
    
    if isnothing(idx)
        throw(ArgumentError("Specified temperature T = $T K not found in the dataset. Available temperatures (approximate): $(round.(unique_temps_all, digits=2))"))
    end
    
    # Get data for the specified temperature
    l_for_T = temp_specific_data_all[idx][1] # loading is the first element
    p_for_T = temp_specific_data_all[idx][2] # pressure is the second element

    # Define the series for plotting
    @series begin
        seriestype := :scatter # Plot data as points
        markersize --> 4
        markerstrokewidth --> 0.5
        
        # Label for the series
        # Since it's always one temperature, a simple "Data" label might suffice,
        # or include the temperature if preferred, especially if combined with model plots.
        label --> Printf.@sprintf("Data at T = %.2f K", unique_temps_all[idx])
        # Or simply: label --> "Data"
            
        # Data for the series: pressure (x), loading (y)
        p_for_T, l_for_T
    end
    # return nothing # Important if all series are generated via @series
end

