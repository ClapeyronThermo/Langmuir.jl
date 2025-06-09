# [Getting Started with Langmuir.jl](@id getting_started)

This is an introductory tutorial for Langmuir.jl. We will demonstrate the basics of the package by building an isotherm model and estimating properties with it.


## Installing Langmuir.jl 

To install Langmuir, use the Julia package manager. This will install the latest release of the package in the Julia package registry.

```julia
using Pkg; Pkg.add("Langmuir")
```

If you want to install the latest development version, you can do so by using the following command:

```julia
using Pkg; Pkg.add(url = "https://github.com/ClapeyronThermo/Langmuir.jl")
```

## Initializing an `IsothermModel` and estimating properties for single component adsorption

In this package, we support several isotherm models. You can refer to the list of supported models here. Here is how you construct a [`LangmuirS1`](@ref) model.

```@example lang1
using Langmuir #hide
M = 1.913 # mol.kg⁻¹
K₀ = 6.82e-10 # Pa⁻¹
E  = -21_976.40 # J.mol⁻¹
isotherm = LangmuirS1(M, K₀, E)
```

You can use a instantiated model to estimate the equilibrium properties of the adsorption system. To estimate the loading (amount of adsorbate per mass of adsorbent) in the adsorbent, given the temperature `T` and pressure `p`, you can do as follows:
 
```@example lang1
p = 101325.0
T = 298.15
l = loading(isotherm, p, T)
```

When estimating loading with a model, it is common to plot isotherms, i.e., pressure vs loading for a fixed temperature. To do it, you can use custom plot call to visualize isotherms using `plot(isotherm, T, P_range)`, where `T` is the temperature you want to analyze, `P_range` is a Tuple specifying minimum and maximum of the range. Below is an example of how to plot two isotherms at two different temperatures.

```@example lang1
using Plots #hide
p_range = (0.0, 1e5)
plot(isotherm, 300.0, p_range)
plot!(isotherm, 400.0, p_range)
```

To finish this section for single component adsorption, one can also estimate the isosteric heat of adsorption by calling `isosteric_heat(model, Vg, p, T)` where Vg is the molar volume of the gas phase, `p` the pressure in Pascal and `T` the temperature in Kelvin. For the Langmuir model, the isosteric heat should be constant and equal to the energy parameter `E`. You can plot the isosteric heat either as a function of the pressure or loading.

 Both the `loading` and `isosteric_heat` functions were made to accept scalar arguments. But they can be used to estimate the loading and isosteric heat for a range of pressures with the broadcasting Julia syntax (you insert a dot after the function name) `loading.(model, p, T)` where `p` is an iterable (Vector or Tuple). 

```@example lang1
P = 0.0:1000.0:1e5 # A range of pressures from 0 to 100 kPa
ΔH = round.(isosteric_heat.(isotherm, P[2:end], 300.0), digits = 8)#Calculates the isosteric heat for a range of pressures
l_at_300 = loading.(isotherm, P[2:end], 300.0) #Calculates the loading at 300 K for a range of pressures since this is a more common visualization
scatter(l_at_300, ΔH, size = (500, 250),  ylabel = "Isosteric heat (J/mol)", xlabel = "loading (mol/kg)", label = "Estimated isosteric heat with AD")
plot!([first(l_at_300), last(l_at_300)], [E, E], label = "Expected value") 
```

## Estimating properties in multicomponent adsorption.

When it comes to estimating properties in multicomponent adsorption, the Ideal Adsorption Solution Theory (IAST) has been proven accurate for a number of systems. It allows one to estimate multicomponent adsorption behavior from single component isotherms.

When formulated, estimating the loading with IAST becomes a nonlinear solve problem which can be solved in different ways. Here, we support the **Nested Loop** and **FastIAS** methods. To know more about the two and which one to choose, refer to this paper: 10.1002/aic.14684.

It can be shown analytically that IAST estimation of multicomponent loading is the same as the extendend Langmuir method when the parameter $M_i$ (saturation loading) are the same for all components, i.e., $M_1 = M_2 = ... = M_{N_c}$. The extended langmuir has the form $n_i = \frac{M_i \times K_{i,0} \exp{\frac{\Delta H}{RT}}}{1 + \sum_i K_i \times P_i}$. Below you can see a numerical verification of IAST for that condition.

```@example multi1
using Langmuir
isotherm_1 = LangmuirS1(1.913, 6.82e-10, -21_976.40)
isotherm_2 = LangmuirS1(1.913, 1.801e-9, -16_925.01)
iastmodel = IASTModels(isotherm_1, isotherm_2)
ext_langmuir = ExtendedLangmuir(isotherm_1, isotherm_2)
loading_1, loading_2 = loading(iastmodel, 101325.0, 300., [0.5, 0.5])


loading_1_expected, loading_2_expected = loading(ext_langmuir, 101325.0, 300., [0.5, 0.5])
 

println("IAST estimated loading for component 1 is: ", round(loading_1, digits = 12))
println("Extended langmuir estimated loading for component 1 is: ", round(loading_1_expected, digits = 12))
```
