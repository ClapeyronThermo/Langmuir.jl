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

Ethane and ethylene data are stored in  `ethane_tpl_data.csv` and `ethylene_tpl.csv` files. To read csv files, many options are available in Julia. Here, we will use `DelimitedFiles.jl`.

In the file, pressure unit is in bar, temperature in Kelvin and loading in mmol/g or mol/kg. 

Let's read and visualize the isotherms at the different temperatures.

```@example fitting
using Plots, DelimitedFiles, Langmuir
ethane_data_path = joinpath(@__DIR__, "sample_data/ethane_tpl_data.csv")
ethylene_data_path = joinpath(@__DIR__, "sample_data/ethylene_tpl_data.csv")
ethane_data = readdlm(ethane_data_path, ',')
P_ethane = ethane_data[:, 2]*1e5
T_ethane = ethane_data[:, 1]
l_ethane = ethane_data[:, 3]
d_ethane = isotherm_data(P_ethane, l_ethane, T_ethane)
Ts_ethane, lp_ethane = split_data_by_temperature(d_ethane)

ethylene_data = readdlm(ethylene_data_path, ',')
P_ethylene = ethylene_data[:, 2]*1e5
T_ethylene = ethylene_data[:, 1]
l_ethylene = ethylene_data[:, 3]
d_ethylene = isotherm_data(P_ethylene, l_ethylene, T_ethylene)
Ts_ethylene, lp_ethylene = split_data_by_temperature(d_ethylene)


scatter(lp_ethane[1][2], lp_ethane[1][1], label = "T = 283.00 K (Ethane)", 
xlabel = "Pressure [Pa]", ylabel = "loading [mol/kg]", m = (4, :white, stroke(1, :slateblue2)),framestyle=:box, markershape = :circle, size = (600, 300))
scatter!(lp_ethane[3][2], lp_ethane[3][1], label = "T = 323.00 K (Ethane)", m = (4, :white, stroke(1, :lightslateblue)), markershape = :circle)

scatter!(lp_ethylene[1][2], lp_ethylene[1][1], label = "T = 283.00 K (Ethylene)", color = :mediumspringgreen, markershape = :square, m = (3, :white, stroke(1, :springgreen2)))
scatter!(lp_ethylene[3][2], lp_ethylene[3][1], label = "T = 323.00 K (Ethylene)",
markershape = :square, m = (3, :white, stroke(1, :springgreen2)))
```

Following the reference manuscript, the quadratic isotherm is the chosen model for fitting the data. In this tutorial, the same strategy is used.

```@example fitting
import Langmuir: x0_guess_fit, to_vec
#Fitting ethane
x0_ethane = to_vec(x0_guess_fit(Quadratic, d_ethane))
lb_ethane = (1e-25, 1e-25, 1e-4, -80_000.0, -80_000.0)
ub_ethane = (1e-1, 1e-1, 100., -1_000.0, -1_000.0)

prob_ethane = IsothermFittingProblem(Quadratic{eltype(d_ethane)}, d_ethane, nothing, abs2, x0_ethane, lb_ethane, ub_ethane) #Bounds have to be manually tweaked. Default interval is too large
alg = DEIsothermFittingSolver(max_steps = 5000, population_size = 50,
logspace = true, verbose = true)
loss_fit_ethane, ethane_isotherm = fit(prob_ethane, alg)
println("Fitting loss for ethane is $loss_fit_ethane")
println(ethane_isotherm)
```

```@example fitting
#Plotting ethane fitting
loading1_ethane = loading_at_T(ethane_isotherm, lp_ethane[1][2], Ts_ethane[1])
loading3_ethane = loading_at_T(ethane_isotherm, lp_ethane[3][2], Ts_ethane[3])
plot!(sort(lp_ethane[1][2]), sort(loading1_ethane), color = :slateblue2, label = "Quadratic - 283.0 K")
plot!(sort(lp_ethane[2][2]), sort(loading3_ethane), color = :slateblue2, label = "Quadratic - 323.0 K")
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
loading1_ethylene = loading_at_T(ethylene_isotherm, lp_ethylene[1][2], Ts_ethylene[1])
loading3_ethylene = loading_at_T(ethylene_isotherm, lp_ethylene[1][2], Ts_ethylene[3])
plot!(sort(lp_ethylene[1][2]), sort(loading1_ethylene), color = :springgreen2, label = "Quadratic - 283.0 K")
plot!(sort(lp_ethylene[1][2]), sort(loading3_ethylene), color = :springgreen2, label = "Quadratic - 323.0 K")
```

From the fitting results, the isotherm parameters closely match those reported in the reference manuscript. For example, the parameter ``M`` was estimated at 3.68 in the manuscript, while it is 3.64 in this analysis. However, the energy parameters exhibit a greater discrepancy. Without calorimetric data, it's challenging to determine the accuracy of these values, as parameter identifiability is limited when relying solely on loading data.

```@example fitting
#Ethane isosteric heat
P_C₂_283K = sort(lp_ethane[1][2][2:end])
l_C₂_283K = sort(loading1_ethane[2:end])
ΔH_C₂_283K = map(p -> isosteric_heat(ethane_isotherm, p, 283.0), P_C₂_283K)

#Ethylene isosteric heat
P_C₂₌_283K = sort(lp_ethylene[1][2][2:end])
l_C₂₌_283K = sort(loading1_ethylene[2:end])
ΔH_C₂₌_283K = map(p -> isosteric_heat(ethylene_isotherm, p, 283.0), P_C₂₌_283K)

scatter(l_C₂_283K, -ΔH_C₂_283K, xlabel = "Loading (mol/kg)", ylabel = "Isosteric Heat (J/mol)", m = (4, :white, stroke(1, :lightslateblue)), markershape = :circle, label = "Ethane - 283.0 K", size = (600, 300))
scatter!(l_C₂₌_283K, -ΔH_C₂₌_283K, xlabel = "Loading (mol/kg)", ylabel = "Isosteric Heat (J/mol)", markershape = :square, m = (3, :white, stroke(1, :springgreen2)), label = "Ethylene - 283.0 K")
```
The plot shows the isosteric heat of adsorption for ethane and ethylene at 283.0 K as a function of loading. Ethane exhibits a steeper decline in adsorption heat, suggesting stronger initial interactions that weaken significantly as loading increases, whereas ethylene's decline is more gradual, indicating a slower reduction in adsorption strength. 


Lastly, the final analysis will utilize IAST to make predictions on the binary adsorption of ethane and ethylene.

```@example fitting
y_C₂ = range(0.0, 1.00, 51) |> collect
success_y_C₂ = []
x_C₂ = []
x_C₂₌ = []
p = 1*101325.0
T = 303.0
for y_C₂_ᵢ in y_C₂
   y = [y_C₂_ᵢ, (1.0 - y_C₂_ᵢ)]
   models = [ethane_isotherm, ethylene_isotherm]
   (n_total, x, is_success) = iast(models, p, T, y, IASTNestedLoop(), maxiters = 1500, abstol = 1e-10, reltol = 1e-10)

   if is_success == :success
      println(is_success)

      push!(x_C₂, x[1])
      push!(x_C₂₌, x[2])
      push!(success_y_C₂, y_C₂_ᵢ)

      else

      nothing

   end
end

plot(success_y_C₂, x_C₂, xlabel = "x", ylabel = "y", label = "ethane", framestyle=:box, size = (600, 300))
plot!(1.0 .- success_y_C₂, x_C₂₌, label = "ethylene")
```

Here you can see the x-y plot for both components.  

