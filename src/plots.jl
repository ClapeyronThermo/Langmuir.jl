"""
# This file contains the plotting recipes for the IsothermModel.

"""

@recipe function f(model::IsothermModel, T::Number, p_range::Tuple{<:Number, <:Number}; npoints = 100)
    # --- Scientific Plot Defaults ---
    # Plot Size and Resolution
    size --> (600, 400) # width, height in pixels
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

    end
