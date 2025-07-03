# [Comparing the estimation of isosteric heat from isotherms with Langmuir.jl](@id isosteric)

During adsorption, heat is released as adsorbate molecules transition to a lower energy state on the adsorbent surface compared to the bulk gas phase. This released heat partially accumulates in the adsorbent, causing a temperature rise on its surface, which in turn can decelerate the adsorption process - as adsorption is an exothermic process.

Predicting equilibrium loading at a given temperature and pressure often receives greater focus. However, accurately predicting the energy release as a function of these same variables is equally critical, as it is as impactful as loading in the adsorption process.

On an energetically heterogeneous surface, the isosteric heat decreases as surface loading increases. Isotherms like the single-site Langmuir model, which assume a constant heat of adsorption regardless of surface loading, are therefore often inadequate for accurately representing experimental data in many cases.

Other isotherms can account for surface heterogeneity such as Toth:

$n = \frac{MKP}{(1 + M(K P)^f)^{1/f}}$ 

However, exhibits an unrealistic infinitely large negative value at high surface loading. below you can see the behavior of the isosteric heat for the toth isotherm as a function of the loading. You can also see the comparison of the analytical expression (which is very tedious to derive) and the one given by `Langmuir.jl` with automatic differentiation.

```@example Toth
using Langmuir, Plots
toth = Toth(7.464, 3.6e-7, 8.3144*-5649.81, 0.5, 50.22)
p_range = 1e-5:500.0:10*101325.0 |> collect
loading_350_t = loading.(toth, p_range, 350.0)
loading_300_t = loading.(toth, p_range, 300.0)
ΔH_350_t = isosteric_heat.(toth, p_range, 350.0)
ΔH_300_t = isosteric_heat.(toth, p_range, 300.0)
plot(loading_300_t/toth.M, ΔH_300_t, framestyle=:box, title = "Toth", xlabel = "covered fraction", ylabel = "ΔH (J/mol)", label = "Automatic Differentiation - 300 K")
plot!(loading_350_t/toth.M, ΔH_350_t, framestyle=:box, label = "Automatic Differentiation - 350 K")

function Q_st1(model::Toth, n1, T)
    E1 = model.E
    β = model.β
    R = 8.31446261815324  # Assuming Rgas gives the gas constant for the type of T
    f = model.f₀ - model.β / T
    n1_0 = model.M  # Assuming saturation loading as reference loading

    # Calculate the two terms inside the brackets in the equation
    term1 = log.(n1 ./ n1_0) ./ (1 .- (n1 ./ n1_0).^f)
    term2 = log.((n1 ./ n1_0) ./ (1 .- (n1 ./ n1_0).^f).^(1/f))

    # Calculate Q_st1 using the main formula
    Q_st1_value = E1 .+ (β * R / f) * (term2 .- term1)
    return Q_st1_value
end

ΔH_analytical = Q_st1(toth, loading_350_t, 350.0)

scatter!(loading_350_t/toth.M, ΔH_analytical,
 label = "Analytical - 350 K", m = (3, :white, stroke(1, :red)))
```

It can be noticed that the isosteric heat for the Toth isotherm blows for high surface coverages. 

Multi-site Langmuir can also account for surface heterogeneity. Below you can see the behavior of the isosteric heat as a function of the surface coverage.

```@example Multisite
using Langmuir, Plots #hide
dualsite = MultiSite(LangmuirS1(2.337, 6.6e-11, 8.3144*-5340.87),
 LangmuirS1(3.490, 3.4e-11, 8.3144*-4273.13))
p_range = 1e-5:500.0:10*101325.0 |> collect
loading_270 = loading.(dualsite, p_range, 270.0)
loading_350 = loading.(dualsite, p_range, 350.0)
ΔH_270 = isosteric_heat.(dualsite, p_range, 270.0)
ΔH_350 = isosteric_heat.(dualsite, p_range, 350.0)
plot(loading_270/(2.337 + 3.490), ΔH_270, framestyle=:box, title = "Dualsite Langmuir", xlabel = "covered fraction", ylabel = "ΔH (J/mol)", label = "270 K")
plot!(loading_350/(2.337 + 3.490), ΔH_350, framestyle=:box, label = "350 K")
```
It can bee seen that the isosteric heat presents an s-shape varying from more energetic to less energetic sites.

[Some literature](https://doi.org/10.1007/s10450-020-00296-3) points that these behaviors are non-physical and potentially problematic when trying to model thermal effects in adsorption.

To overcome it, the adsorption Nonrandom Two-Liquid (aNRTL) activity coefficient model into an activity-based formulation for Langmuir isotherm, [Chang et al. 2020](https://doi.org/10.1007/s10450-019-00185-4) proposed a thermodynamic Langmuir (tL) which seems to have superior properties for predicting the isosteric heat compared to Toth and Multisite Langmuir.

The equation for the loading $n_i$ in terms of $\gamma_i$ and $\gamma_\phi$ is:

$n_i = \frac{M K_i P}{\frac{\gamma_i}{\gamma_\phi} + K_i P}$

where $K_i$ is also a function of temperature $K_i = K_i^\circ\exp{\frac{-E}{RT}}$, and the activity coefficients $\gamma_i$ and $\gamma_\phi$ are given by:

$\gamma_i = \exp\left(\frac{d \frac{g^E}{RT}}{d \theta_i}\right)$

$\gamma_\phi = \exp\left(\frac{d \frac{g^E}{RT}}{d \theta_{\phi}}\right)$

The Gibbs excess free energy term, $\frac{g^E}{RT}$, is expressed as:

$\frac{g^E}{RT} = \frac{\theta_i \theta_\phi \tau_{i\phi} (G_{i\phi} - 1)}{\theta_i G_{i\phi} + \theta_\phi}$

where, $B_{i\phi}$ is a model parameter, $\tau_{i\phi} = B_{i\phi} / T$, $G_{i\phi} = \exp(-0.3 \cdot \tau_{i\phi})$.

Here, $\theta_i$ and $\theta_\phi$ are the coverage terms for the adsorbed species and adsorption sites, 
respectively ($\theta_i  + \theta_{\phi} = 1$), and $T$ is the temperature.

Below you can see how to iniatilize the thermodynamic langmuir model in `Langmuir.jl` and use it to predict the loading and isosteric heat.

```@example tLangmuir
using Langmuir, Plots #hide
tlang = ThermodynamicLangmuir(5.890, 6.1e-11, -4599.86*8.3144, -762.51)
p_range = 1e-5:500.0:10*101325.0 |> collect
loading_270 = loading.(tlang, p_range, 270.0)
loading_350 = loading.(tlang, p_range, 350.0)
ΔH_270 = isosteric_heat.(tlang, p_range, 270.0)
ΔH_350 = isosteric_heat.(tlang, p_range, 350.0)
plot(loading_270/tlang.M, ΔH_270, framestyle=:box, xlabel = "covered fraction", ylabel = "ΔH (J/mol)", title = "tLangmuir", label = "270 K")
plot!(loading_350/tlang.M, ΔH_350, framestyle=:box, label = "350 K")
plot!([minimum(loading_270), maximum(loading_270)]./tlang.M,  [-4599.86*8.3144, -4599.86*8.3144], label = "Langmuir")
```

It is evident that the thermodynamic Langmuir model exhibits behavior that is notably different from both the Toth and multisite Langmuir models. [As highlighted in existing literature](https://link.springer.com/article/10.1007/s10450-019-00185-4), the thermodynamic Langmuir model offers more accurate estimates of isosteric heat while maintaining high predictive accuracy for equilibrium loading. It can also be use to predict binary adsorption. 

Analytical expressions for the isosteric heat can become quite extensive. In fact, there is an [entire manuscript](https://aiche.onlinelibrary.wiley.com/doi/10.1002/aic.17186) dedicated to deriving such expressions for a number of isotherms. 

In `Langmuir.jl`, we leverage automatic differentiation to get accurate estimates of the required derivatives for the isosteric heat without requiring such extensive derivations.

