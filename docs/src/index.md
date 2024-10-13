```@meta
CurrentModule = Langmuir
```

# Langmuir.jl 

Langmuir.jl is a powerful Julia library designed to model adsorption equilibria for both single and multi-component systems.

For single-component adsorption, the library offers a wide range of isotherms, from simple one-parameter models like Henry's law to more complex two- and three-parameter models such as the Langmuir (Single and MultiSite), Freundlich, Temkin, Redlich-Peterson, Toth, and Sips isotherms. These models include temperature-dependent parameters, which are essential for estimating the isosteric heat of adsorption from pressure-loading data sets at varying temperatures. For a complete list of available models, refer to ___.

For multicomponent systems, Langmuir.jl employs the Ideal Adsorption Solution Theory (IAST) to predict the loading of different components, given their bulk pressure, temperature and individual isotherms. This approach allows for predictive modeling of multicomponent adsorption.

A complete workflow for using Langmuir.jl typically includes the following steps:

   - **Data Preparation:** Load your dataset, which should include adsorption loading, pressure, temperature, and heat data.

   - **Isotherm Fitting:** Fit the appropriate isotherm model to your data using a global optimization method to ensure accurate parameter estimation.

   - **Multicomponent Loading** Estimation: Use Ideal Adsorption Solution Theory (IAST) to predict the loading of each component in a multicomponent system based on bulk pressure and temperature.

   - **Isosteric Heat Calculation:** Estimate the isosteric heat of adsorption for single or multicomponent systems using the library's built-in functions, which account for temperature dependencies in the isotherm parameters.


### Authors

- [Andrés Riedemann](mailto:andres.riedemann@gmail.com), University of Concepción
- [Vinicius Santana](mailto:vinicius.viena1@gmail.com), Norwegian University of Science and Technology - NTNU.
- Pierre Walker - Caltech. 
