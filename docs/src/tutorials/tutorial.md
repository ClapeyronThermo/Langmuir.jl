# [A typical workflow with Langmuir.jl](@id Tutorial)

In this tutorial, we will go through a typical workflow for processing and analyzing adsorption equilibrium data:

1. **Loading adsorption data** for the pure components in the mixture of interest.
2. **Fitting an isotherm** for each component.
3. **Analyzing the fitting quality** by:
   - Plotting predictions vs. experimental data.
   - Performing residual analysis.
4. **Estimating single-component properties** such as the isosteric heat of adsorption.
5. **Estimating multicomponent adsorption loading** using the Ideal Adsorbed Solution Theory (IAST).


For this tutorial, we will work with a binary system consisting of ethane and ethylene. The goal is to investigate their adsorption properties in DUT-8, a material known for being ethane-selective.

This case study is based on the research presented in the following paper:

- **Reference:** Santana, V. V., Carmo, P., Seabra, R., Rodrigues, A. E., Ribeiro, A. M., Nogueira, I. B. R., Yoon, J. W., Cho, K. H., Yun, J. S., Lee, U.-H., Kim, K., Ferreira, A. F. P. (2024). *Ethylene Purification by Pressure Swing Adsorption with the Paraffin Selective Metal–Organic Framework─DUT-8*. **Industrial & Engineering Chemistry Research**, 63(5), 2307–2319. [DOI: 10.1021/acs.iecr.3c02808](https://doi.org/10.1021/acs.iecr.3c02808).

Ethane and ethylene data are stored in  `ethane_tpl_data.csv` and `ethylene_tpl.csv` files. To read csv files, many options are available in Julia. Here, we will use `DelimitedFiles.jl`. We also use the `isotherm_data` function to wrap the data into a format suitable for processing the data with Langmuir.jl. The data files contain pressure, loading, and temperature information for each component. **Note that keeping this order is crucial** for the `isotherm_data` function to work correctly.

In the file, pressure unit is in bar, temperature in Kelvin and loading in mmol/g or mol/kg. 

Let's read and visualize the isotherms at the different temperatures. You can use a custom plotting function to visualize the data by calling `plot(data, temperature)`. This function will plot the pressure against the loading for a given temperature.

```@example fitting
using Plots, DelimitedFiles, Langmuir
ethane_data_path = joinpath(@__DIR__, "sample_data/ethane_tpl_data.csv")
ethylene_data_path = joinpath(@__DIR__, "sample_data/ethylene_tpl_data.csv")
ethane_data = readdlm(ethane_data_path, ',')
P_ethane = ethane_data[:, 2]*1e5 # Convert bar to Pa
T_ethane = ethane_data[:, 1]
l_ethane = ethane_data[:, 3]
d_ethane = isotherm_data(P_ethane, l_ethane, T_ethane)

ethylene_data = readdlm(ethylene_data_path, ',')
P_ethylene = ethylene_data[:, 2]*1e5 # Convert bar to Pa
T_ethylene = ethylene_data[:, 1]
l_ethylene = ethylene_data[:, 3]
d_ethylene = isotherm_data(P_ethylene, l_ethylene, T_ethylene) #Alwas read in order of Pressure, Loading, Temperature

figsize = (500, 500/1.618)
plot(d_ethane, 283.0, label = "Ethane at 283K", m = (3, :white, stroke(1, :blue)), size = figsize, xlabel = "Pressure (Pa)", ylabel = "Loading (mol/kg)", markershape = :circle)
plot!(d_ethane, 323.0, label = "Ethane at 323K", markershape = :square, m = (3, :white, stroke(1, :blue)))

plot!(d_ethylene, 283.0, label = "Ethylene at 283K", markershape = :circle, 
m = (3, :white, stroke(1, :red)))
scatter!(d_ethylene, 323.0, label = "Ethylene at 323K",
markershape = :square, m = (3, :white, stroke(1, :red)))
```

Following the reference manuscript, the quadratic isotherm is the chosen model for fitting the data. In this tutorial, the same strategy is used. Note that the bounds for the parameters were be manually set since the default ones were too large. Once the bounds are custom, you also need to provide the initial guess. We have default initial guesses for each isotherm, but to expose it, you need to `import` the `x0_guess_fit` function. We also import the `to_vec` function to convert the initial guess into a vector format, which is required by the fitting problem.

 Also observe that the argument next to d_ethane is `nothing`. This is because the fitting problem does not involve any additional data, such as calorimetric data, to estimate the isotherm parameters.

```@example fitting
import Langmuir: x0_guess_fit, to_vec
#Fitting ethane
x0_ethane = to_vec(x0_guess_fit(Quadratic, d_ethane))
lb_ethane = (1e-25, 1e-25, 1e-4, -80_000.0, -80_000.0)
ub_ethane = (1e-1, 1e-1, 100., -1_000.0, -1_000.0)

prob_ethane = IsothermFittingProblem(Quadratic{eltype(d_ethane)}, d_ethane, nothing, abs2, x0_ethane, lb_ethane, ub_ethane) #Bounds have to be manually tweaked. Default interval is too large
alg = DEIsothermFittingSolver(max_steps = 500, population_size = 150,
logspace = true, verbose = true)
loss_fit_ethane, ethane_isotherm = fit(prob_ethane, alg)
println("Fitting loss for ethane is $loss_fit_ethane")
println(ethane_isotherm)
```

The fitting results can be visualized by plotting the isotherm predictions against the experimental data. The `plot!` function can be used to overlay the fitted isotherm on the experimental data. This is also a custom case of the `plot` function from the `Plots.jl` to facilitate the visualization of isotherms. You can use it to plot the isotherm predictions for a specific temperature and pressure range. 


```@example fitting
#Plotting ethane fitting
plot!(ethane_isotherm, 283.0, (0.0, maximum(P_ethane)), 
color = :blue, label = "Quadratic Ethane - 283.0 K", linestyle = :dash)
plot!(ethane_isotherm, 323.0, (0.0, maximum(P_ethane)), 
color = :blue, label = "Quadratic Ethane - 323.0 K")
```

```@example fitting
#Fitting Ethylene
x0_ethylene = to_vec(x0_guess_fit(Quadratic, d_ethylene))
lb_ethylene = (1e-25, 1e-25, 1e-4, -80_000.0, -80_000.0)
ub_ethylene = (1e-1, 1e-1, 100., -500.0, -500.0)
prob_ethylene = IsothermFittingProblem(Quadratic{eltype(d_ethylene)}, d_ethylene, nothing, abs2, x0_ethylene, lb_ethylene, ub_ethylene)
loss_fit_ethylene, ethylene_isotherm = fit(prob_ethylene, alg)
println("Fitting loss for ethane is $loss_fit_ethylene")
println(ethylene_isotherm)
```

```@example fitting
#Plotting ethylene isotherms
plot!(ethylene_isotherm, 283.0, (0.0, maximum(P_ethane)),
label = "Quadratic Ethylene - 283.0 K", linestyle = :dash, color = :red)
plot!(ethylene_isotherm, 323.0, (0.0, maximum(P_ethane)),
label = "Quadratic Ethylene - 283.0 K", color = :red)
```

From the fitting results, the isotherm parameters closely match those reported in the reference manuscript. For example, the parameter `M` was estimated at 3.68 in the manuscript, while it is 3.64 in this analysis. However, the energy parameters exhibit a greater discrepancy. Without calorimetric data, it's challenging to determine the accuracy of these values, as parameter identifiability is limited when relying solely on loading data.

```@example fitting
#Ethane isosteric heat
is_283K = d_ethane.T .== 283.0
P_C₂_283K = d_ethane.p[is_283K]
l_C₂_283K = d_ethane.l[is_283K]
ΔH_C₂_283K = isosteric_heat.(ethane_isotherm, P_C₂_283K, 283.0)

#Ethylene isosteric heat
is_283K_ethylene = d_ethylene.T .== 283.0
P_C₂₌_283K = d_ethylene.p[is_283K_ethylene]
l_C₂₌_283K = d_ethylene.l[is_283K_ethylene]
ΔH_C₂₌_283K = isosteric_heat.(ethylene_isotherm, P_C₂₌_283K, 283.0) 

scatter(l_C₂_283K, ΔH_C₂_283K, xlabel = "Loading (mol/kg)", ylabel = "Isosteric Heat (J/mol)", m = (4, :white, stroke(1, :lightslateblue)), markershape = :circle, label = "Ethane - 283.0 K", size = (600, 300))
scatter!(l_C₂₌_283K, ΔH_C₂₌_283K, xlabel = "Loading (mol/kg)", ylabel = "Isosteric Heat (J/mol)", markershape = :square, m = (3, :white, stroke(1, :springgreen2)), label = "Ethylene - 283.0 K")
```
The plot shows the isosteric heat of adsorption for ethane and ethylene at 283.0 K as a function of loading. Ethane exhibits a steeper decline in adsorption heat, suggesting stronger initial interactions that weaken significantly as loading increases, whereas ethylene's decline is more gradual, indicating a slower reduction in adsorption strength. 

Lastly, the final analysis will utilize IAST to make predictions on the binary adsorption of ethane and ethylene. There are two ways to use IAST in Langmuir.jl. 
You can directly use the `iast` function, which gives you access to extra solution parameters such as convergence criteria, total loading and composition. Another way is to wrap the isotherms into the `IASTModels` type and use it to estimate quantities such as the loading and the isosteric heat. Below, we are using the first approach, which is more straightforward for this case.

```@example fitting
y_C₂ = range(0.0, 1.00, 51) #51 points from 0.0 to 1.0
success_y_C₂ = []
x_C₂ = []
x_C₂₌ = []
p = 1*101325.0
T = 303.0
for y_C₂_ᵢ in y_C₂
   y = [y_C₂_ᵢ, (1.0 - y_C₂_ᵢ)]
   models = [ethane_isotherm, ethylene_isotherm]
   (n_total, x, is_success) = iast(models, p, T, y, FastIAS(), maxiters = 1500, abstol = 1e-10, reltol = 1e-10)

   if is_success == :success
      println(is_success)

      push!(x_C₂, x[1])
      push!(x_C₂₌, x[2])
      push!(success_y_C₂, y_C₂_ᵢ)

      else

      nothing

   end
end

plot(x_C₂, success_y_C₂, xlabel = "x (adsorbed phase)", ylabel = "y (gas phase)", label = "ethane", framestyle=:box, size = (600, 300))
plot!(x_C₂₌, 1.0 .- success_y_C₂, label = "ethylene")
```
 

