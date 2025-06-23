# [Finding the best isotherm for your data with Langmuir.jl](@id best isotherm)

Working with single component models is an important part of adsorption modeling, and it is often difficult to choose the best model for a given system. 

To illustrate the investigation of choosing the right isotherm using `Langmuir.jl`, this tutorial was prepared. The data was obtained from [Möllmer et al. (2012)](@cite Mollmer2012) where isotherm data of methane adsorption on
Cu(Me-4py-trz-ia) material is discussed. The data was extracted using an automated symbol extraction tool from the original publication and may containg some errors, specially at low pressures where reading is difficult as the markers overlap.

The processing of the data consists of the following steps:
1. **Read the data**: The data is read from CSV files, which contain pressure and loading values for each component at different temperatures.
2. **Create isotherm data structures**: The data is converted into `IsothermData` structures, which are used to store the pressure, loading, and temperature information.
3. **Merge isotherm data**: The individual isotherm data for each temperature is merged into a single structure using the `merge_isotherm_data` function.
4. **Plot the isotherm data**: The isotherm data can then be plotted to visualize the adsorption behavior of each component at different temperatures.

```@example multicomponent
using DelimitedFiles, Plots, Langmuir
Temperatures = [273.0, 298.0, 323.0]

#Data only has pressure and loading, so we need to create a temperature vector
CH4_273 = readdlm(joinpath(@__DIR__, "sample_data/ch4_273K.csv"), ',')
T⃗_237 = fill(273.0, size(CH4_273, 1))

CH4_298 = readdlm(joinpath(@__DIR__, "sample_data/ch4_298K.csv"), ',')
T⃗_298 = fill(298.0, size(CH4_298, 1))

CH4_323 = readdlm(joinpath(@__DIR__, "sample_data/ch4_323K.csv"), ',')
T⃗_323 = fill(323.0, size(CH4_323, 1))

#Individual isotherm data for each temperature
CH4_data_273 = isotherm_data(CH4_273[:, 1], CH4_273[:, 2], T⃗_237) # Pressure, Loading, Temperature
CH4_data_298 = isotherm_data(CH4_298[:, 1], CH4_298[:, 2], T⃗_298)
CH4_data_323 = isotherm_data(CH4_323[:, 1], CH4_323[:, 2], T⃗_323)

# Merge the isotherm data into a single structure
CH4_data = merge_isotherm_data(CH4_data_273, CH4_data_298, CH4_data_323)

fig1 = Plots.plot()
xlabel = "Pressure (MPa)"
ylabel = "Loading (mol/kg)"
plot!(fig1, CH4_data, 273.0, label = "CH4 at 273K", markershape = :hexagon, xlabel = xlabel, ylabel = ylabel, m = (4, :white, stroke(1, :blue)), size = (800, 350), legend_columns=1)
plot!(fig1, CH4_data, 298.0, label = "CH4 at 298K", markershape = :square, m = (4, :white, stroke(1, :green)))
plot!(fig1, CH4_data, 323.0, label = "CH4 at 323K", markershape = :circle, m = (4, :white, stroke(1, :red)))
```

!!! note
    Observe that the data set automatically includes a variance column, which refers to the variance in the loading measuraments and is set to 1.0 for all data points. This is a placeholder and can be replaced with actual variance values if available. If you are less certain about some measurements than others, you can set the variance to a higher value for those points, which will make the fitting algorithm less sensitive to them.

The reference manuscript fitted a Toth model to CH4. Let's try a few models to see which one fits the data best. The `Langmuir` package provides a variety of isotherm models, including `LangmuirS1`, `Toth`, `Sips` and `Thermodynamic Langmuir`. 

Fitting a Langmuir model is quite simple as the parameters are easy to initialize. We will use it as an initial guess for all models, but we want to compare first the Langmuir single-site model with the Thermodynamic Langmuir model, which is an implicit model that takes into consideration an activity coefficient.  

```@example multicomponent
using Langmuir
alg = DEIsothermFittingSolver(max_steps = 500, population_size = 100,
logspace = true, verbose = true)
prob_lang = IsothermFittingProblem(LangmuirS1, CH4_data, abs2)
loss_lang, model_lang = fit(prob_lang, alg)
plot!(fig1, model_lang, 273.0, (0.0, maximum(CH4_data.p)), color = :blue)
plot!(fig1, model_lang, 298.0, (0.0, maximum(CH4_data.p)), color = :green)
plot!(fig1, model_lang, 323.0, (0.0, maximum(CH4_data.p)), color = :red)
```

Now, we will compare similar models, i.e., the Toth and Sips models, which are both more complex than the Langmuir model. The Toth model is a more flexible model that can fit a wider range of isotherm shapes, while the Sips model is a combination of the Langmuir and Freundlich models.

```@example multicomponent 
x0_tlang = [model_lang.M, model_lang.K₀, model_lang.E, -1e-10]
lb_tlang = (1e-10, 1e-25, -70_000.0, -2000.0)
ub_tlang = (20.0, 1e-1, -1000.0, 0.0)
calorimetric_data = nothing # No calorimetric data available for this example
prob_tlang = IsothermFittingProblem(ThermodynamicLangmuir, CH4_data, calorimetric_data, abs2, x0_tlang, lb_tlang, ub_tlang)
alg = DEIsothermFittingSolver(max_steps = 1000, population_size = 500,
logspace = true, verbose = true, time_limit = 15)
loss_tlang, model_tlang = fit(prob_tlang, alg)
plot!(fig1, model_tlang, 273.0, (0.0, maximum(CH4_data.p)), color = :blue, linestyle = :dot)
plot!(fig1, model_tlang, 298.0, (0.0, maximum(CH4_data.p)), color = :green, linestyle = :dot)
plot!(fig1, model_tlang, 323.0, (0.0, maximum(CH4_data.p)), color = :red, linestyle = :dot)
```

The Toth model is more trick to fit as the `f` parameter is made temperature dependent. The default bounds for the parameters are too large, so we need to set them manually. The fitted Langmuir above can be used as a good initial guess for the parameters.

!!! note
    **It is important to set a time limit for the optimization algorithm** to prevent it from running indefinitely in search of a better fit - from previous experience, you won't need more than 10s to fit an isotherm with a good initial guess. You can do this by setting the `time_limit` parameter in the `DEIsothermFittingSolver` constructor.

```@example multicomponent
fig2 = Plots.plot()
xlabel = "Pressure (MPa)"
ylabel = "Loading (mol/kg)"
plot!(fig2, CH4_data, 273.0, label = "CH4 at 273K", markershape = :hexagon, xlabel = xlabel, ylabel = ylabel, m = (4, :white, stroke(1, :blue)), size = (800, 350), legend_columns=1)
plot!(fig2, CH4_data, 298.0, label = "CH4 at 298K", markershape = :square, m = (4, :white, stroke(1, :green)))
plot!(fig2, CH4_data, 323.0, label = "CH4 at 323K", markershape = :circle, m = (4, :white, stroke(1, :red)))


x0_toth = [model_lang.M, model_lang.K₀, model_lang.E, 1.0, -1e-10]
lb_toth = (1e-10, 1e-25, -70_000.0, 0.0, -300.0)
ub_toth = (20.0, 1e-1, -1000.0, 2.0, 0.0)
calorimetric_data = nothing # No calorimetric data available for this example
prob_toth = IsothermFittingProblem(Toth, CH4_data, calorimetric_data, abs2, x0_toth, lb_toth, ub_toth)
alg = DEIsothermFittingSolver(max_steps = 500, population_size = 100,
logspace = true, verbose = true, time_limit = 15)
loss_toth, model_toth = fit(prob_toth, alg)
plot!(fig2, model_toth, 273.0, (0.0, maximum(CH4_data.p)), color = :blue, linestyle = :dash)
plot!(fig2, model_toth, 298.0, (0.0, maximum(CH4_data.p)), color = :green, linestyle = :dash)
plot!(fig2, model_toth, 323.0, (0.0, maximum(CH4_data.p)), color = :red, linestyle = :dash)
```

A good reminder about our formulation of the Toth model is that the $f$ parameter has a temperature dependency of $f = f_0 - \frac{\beta}{T}$. The $\beta$ parameter is initially defined with bounds in $(-\infty, \infty)$ while $f_0$ is bounded to $(0, \infty)$. However, it is good practice to investigate the direction in which $\beta$ value changes with temperature. After empirical testing, it became evident that $\beta$ should be constrained to the interval $(-\infty, 0)$ instead of including positive numbers and the bounds were updated accordingly. Nevertheless, the fitted model tends to push $\beta$ to the lower bound, effectively driving the $f_0$ term towards zero, i.e., the temperature dependency of $f$ that leads to the lowest fitting error effectively becomes $f = -\frac{\beta}{T}$. 

This behavior reflects a common issue in isotherm models with a high number of parameters: parameter correlation leads to identifiability problems, where different parameter combinations can yield similar loss values. As a result, optimization may favor extreme values that are not physically meaningful. 

```@example multicomponent
x0_sips = [model_lang.M, model_lang.K₀, model_lang.E, 1.0, -1.0]
lb_sips = (1e-10, 1e-25, -70_000.0, 0.0, -300.0)
ub_sips = (20.0, 1e-1, -1000.0, 2.0, 0.0)
prob_sips = IsothermFittingProblem(Sips, CH4_data, calorimetric_data, abs2, x0_sips, lb_sips, ub_sips)
alg = DEIsothermFittingSolver(max_steps = 500, population_size = 100,
logspace = true, verbose = true, time_limit = 15)
loss_sips, model_sips = fit(prob_sips, alg)
plot!(fig2, model_sips, 273.0, (0.0, maximum(CH4_data.p)), color = :blue, linestyle = :dashdot)
plot!(fig2, model_sips, 298.0, (0.0, maximum(CH4_data.p)), color = :green, linestyle = :dashdot)
plot!(fig2, model_sips, 323.0, (0.0, maximum(CH4_data.p)), color = :red, linestyle = :dashdot)
```

You can build a table with the results of the fitting to compare the models. The `parameters` field of the `IsothermModel` structure contains the fitted parameters, and the `loss` variable contains the loss of the fitting, i.e., $\mathcal{L} = \sum_i^N \frac{(q_i - q*_i)^2}{\sigma^2_i}$. 

# Langmuir models (Single Site and Thermodynamic)

```@example multicomponent
using DataFrames, Printf
langmuir_results = DataFrame(
    Model = ["Langmuir", "Thermodynamic Langmuir"],
    Loss = round.([loss_lang, loss_tlang], digits=4),
    M = round.([model_lang.M, model_tlang.M], digits=3),
    K₀ = [model_lang.K₀, model_tlang.K₀],
    E = round.([model_lang.E, model_tlang.E], digits=1),
    Bᵢᵩ = [0.0, round(model_tlang.Bᵢᵩ, sigdigits=3)]
)
```

# Complex models table (Toth and Sips)

```@example multicomponent
complex_results = DataFrame(
    Model = ["Toth", "Sips"],
    Loss = round.([loss_toth, loss_sips], digits=4),
    M = round.([model_toth.M, model_sips.M], digits=3),
    K₀ = [model_toth.K₀, model_sips.K₀],
    E = round.([model_toth.E, model_sips.E], digits=1),
    f₀ = [round(model_toth.f₀, digits=3), round(model_sips.f₀, digits=3)],
    β = [round(model_toth.β, sigdigits=3), round(model_sips.β, sigdigits=3)]
);
```

As can be seen, the Toth model achieves the lowest mean squared error. However, it can be seen that all measurements are assigned the same variance, which is not the case in reality. Low pressure measurements are more uncertain than high pressure ones, so we can assign a higher variance to the low pressure points and see the consequences in the fitting process. You can do this by modifying the `σ²` field of the `IsothermData` structure.

