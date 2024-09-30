# [Getting Started with AdsorbedSolutionTheory.jl](@id getting_started)

This is an introductory tutorial for AdsorbedSolutionTheory.jl (AST). We will demonstrate the basics of the package by building an isotherm model and estimating properties with it.


## Installing AdsorbedSolutionTheory.jl 

To install AdsorbedSolutionTheory, use the Julia package manager.

```julia
using Pkg; Pkg.add("AdsorbedSolutionTheory")
```

## Initializing an `IsothermModel` and estimating properties for single component adsorption

In this package, we support several isotherm models. You can refer to the list of supported models here. Here is how you construct a [`Langmuir`](@ref) model.

```@example lang1
using AdsorbedSolutionTheory #hide
M = 1.913 # mol.kg⁻¹
K₀ = 6.82e-10 # Pa⁻¹
E  = -21_976.40 # J.mol⁻¹
isotherm = Langmuir(M, K₀, E)
```

You can use a instantiated model to estimate the equilibrium properties of the adsorption system. To estimate the loading (amount of adsorbate per mass of adsorbent) in the adsorbent, given the temperature `T` and pressure `p`, you can do as follows:
 
```@example lang1
p = 101325.0
T = 298.15
l = loading(isotherm, p, T)
```

When estimating loading with a model, it is common to plot isotherms, i.e., pressure vs loading for a fixed temperature. To do it, you can use the `loading_at_T(isotherm, P, T)` function.

```@example lang1
using Plots #hide
P = 0.0:5_000.0:100_000.0 |> collect
l_at_300 = loading_at_T(isotherm, P, 300.)
l_at_350 = loading_at_T(isotherm, P, 350.)
plot(P, l_at_300, size = (500, 250), label = "300K")
plot!(P, l_at_350, label = "350K")
xlabel!("P (Pa)")
ylabel!("l (mol/kg)")
```

You can also estimate other properties from the isotherm such as the henry coefficient at a given temperature by calling `henry_coefficient(model::IsothermModel, T)`. The henry coefficient should correspond to the slope of the isotherm when $P \rightarrow 0.0$. In AdsorbedSolutionTheory.jl, this is obtained using automatic differentiation and introduces no numerical error in the estimate. You can see in the example below how to visualize the tangent line built from the henry coefficient at $300K$.

```@example lang1
P_ = P[1:3]
plot(P_, l_at_300[1:3], size = (500, 250), label = "300K") 
H = henry_coefficient(isotherm, 300.0)
plot!(P_, H*P_, label = "Tangent line")
```

To finish this section for single component adsorption, one can also estimate the isosteric heat of adsorption by calling `isosteric_heat(model, Vg, p, T)` where Vg is the molar volume of the gas phase, `p` the pressure in Pascal and `T` the temperature in Kelvin. For the Langmuir model, the isosteric heat should be constant and equal to the energy parameter `E`. You can plot the isosteric heat either as a function of the pressure or loading.

Below it is assumed that the ideal gas law is a good approximation to describe the molar volume of the gas phase.

```@example lang1
import AdsorbedSolutionTheory: Rgas
Vg = Rgas(isotherm)*300.0./P[2:end]
ΔH = map(Vg_P -> isosteric_heat(isotherm, first(Vg_P), last(Vg_P), 300.), zip(Vg, P[2:end])) |> x -> round.(x, digits = 7)
scatter(l_at_300[2:end], ΔH, size = (500, 250),  ylabel = "Isosteric heat (J/mol)", xlabel = "loading (mol/kg)", label = "Estimated isosteric heat with AD")
plot!([first(l_at_300), last(l_at_300)], [-E, -E], label = "Expected value") 
```

## Estimating properties in multicomponent adsorption.

When it comes to estimating properties in multicomponent adsorption, the Ideal Adsorption Solution Theory (IAST) has been proven accurate for a number of systems. It allows one to estimate multicomponent adsorption behavior from single component isotherms.

When formulated, estimating the loading with IAST becomes a nonlinear solve problem which can be solved in different ways. Here, we support the **Nested Loop** and **FastIAS** methods. To know more about the two and which one to choose, refer to this paper: 10.1002/aic.14684.

It can be shown analytically that IAST estimation of multicomponent loading is the same as the extendend Langmuir method when the parameter $M_i$ (saturation loading) are the same for all components, i.e., $n_i = \frac{M_i \times K_{i,0} \exp{\frac{\Delta H}{RT}}}{1 + \sum_i K_i \times P_i}$. Below you can see a numerical example of it.

```@example multi1
using AdsorbedSolutionTheory
import AdsorbedSolutionTheory: Rgas
isotherm_1 = Langmuir(1.913, 6.82e-10, -21_976.40)
isotherm_2 = Langmuir(1.913, 1.801e-9, -16_925.01)
models = (isotherm_1, isotherm_2)
(n_total, x, is_success) = iast(models, 101325.0, 300., [0.5, 0.5], FastIAS())
loading_1 =  n_total*x[1]
loading_2 = n_total*x[2]

K_1 = 6.82e-10*exp(21_976.40/Rgas(isotherm_1)/300.)
K_2 = 1.801e-9*exp(16_925.01/Rgas(isotherm_2)/300.)
p_1 = 0.5*101325.0
p_2 = 0.5*101325.0
loading_1_expected = 1.913*K_1*p_1/(1.0 + K_1*p_1 + K_2*p_2)
loading_2_expected = 1.913*K_2*p_2/(1.0 + K_1*p_1 + K_2*p_2)

println("IAST estimated loading for component 1 is: ", round(loading_1, digits = 4))
println("Extende langmuir estimated loading for component 1 is: ", round(loading_1_expected, digits = 4))
```
