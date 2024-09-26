# [Getting Started with AdsorbedSolutionTheory.jl](@id getting_started)

This is an introductory tutorial for AdsorbedSolutionTheory.jl (AST). We will demonstrate the basics of the package by building an isotherm model and estimating properties with it.


## Installing AdsorbedSolutionTheory.jl 

To install AdsorbedSolutionTheory, use the Julia package manager.

```julia
using Pkg; Pkg.add("AdsorbedSolutionTheory")
```

## Initializing an `IsothermModel` and estimating properties

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

When estimating loading with a model, it is common to plot isotherms, i.e., pressure vs loading
for a fixed temperature. To do it, you can use the `loading_at_T` function.

```@example lang1
using Plots #hide
P = 10_000:5_000:100_000
l_at_300 = loading_at_T(isotherm, P, 300.)
l_at_350 = loading_at_T(isotherm, P, 350.)
plot(collect(P), l_at_300, size = (500, 250), label = "300K")
plot!(collect(P), l_at_350, label = "350K")
xlabel!("P (Pa)")
ylabel!("l (mol/kg)")
```


