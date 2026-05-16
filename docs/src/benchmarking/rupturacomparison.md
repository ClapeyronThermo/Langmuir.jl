# Benchmarking

Adsorption equilibrium calculations have been a field of active research and development for decades, and there are many well-implemented open-source libraries that have helped students, researchers, and industrial practitioners achieve their goals in this area for a long time. A few notable examples include PyIAST (Python), Ruptura (C++ with Python bindings), and PyGAPS (Python, which also calls PyIAST for IAST calculations). There is therefore no reason to claim that Langmuir.jl is the only or best solution available, but it fills its own space and offers its own set of features within the ecosystem.

That said, despite the adsorption community being well served by existing tools, Langmuir.jl was born from the belief that, by leveraging Julia and its combination of high performance and friendly syntax, practitioners can expand their toolset for adsorption thermodynamics calculations and tackle problems that were previously out of reach.

With this in mind, this section has two goals. The first is to demonstrate that Langmuir.jl is a reliable library by benchmarking its most important calculations against results from credible sources, namely Ruptura and PyIAST. The second is to show that its performance is on par with Ruptura, which we consider to represent the state of the art in terms of computational performance. Ruptura is a software library for adsorption thermodynamics calculations and breakthrough curve simulation written in C++ with Python bindings, while PyIAST is a Python-native library for adsorption thermodynamics calculations.

## Example 1 - Equimolar 4-component xylene mixture prediction with IAST

The example was taken from the [Ruptura repository](https://github.com/iRASPA/RUPTURA/blob/main/examples/CoBDP-xylenes/iast/run.py) an equimolar mixture of 4 xylenes – o-xylene, m-xylene, p-xylene and ethylbenzene ($25\%$ molar). Pure components isotherms are modeled with dual-site Langmuir-Freundlich isotherms:

$$q_i = \sum_{i = 1}^2 \frac{q^{\textrm{sat}}_i \times K_i \times p^f_i}{1 + K_i \times p^f_i}$$

Ruptura does not account for temperature dependency of parameters as Langmuir.jl. It's not clear why it is the case, but it is also the way that PyIAST operates with the implementation of the isotherms. The parameters at $443.0 K$ are:

| Component    | y    | q_sat,1 | b₁        | n₁   | q_sat,2 | b₂        | n₂  |
|--------------|------|---------|-----------|------|---------|-----------|-----|
| o-xylene     | 0.25 | 4.8     | 1.28×10⁻⁵ | 1.4  | 1.5     | 1.6×10⁻⁴  | 0.7 |
| m-xylene     | 0.25 | 4.8     | 3.19×10⁻⁵ | 1.27 | 1.2     | 4.72×10⁻⁵ | 0.7 |
| p-xylene     | 0.25 | 4.5     | 2.3×10⁻⁶  | 1.7  | 1.6     | 1.46×10⁻⁴ | 0.7 |
| ethylbenzene | 0.25 | 4.5     | 7.38×10⁻⁶ | 1.5  | 1.4     | 1.6×10⁻⁴  | 0.7 |

### Ruptura Syntax

To show how the syntax of Ruptura compares to Langmuir.jl, we display here a code snippet taken from the documented examples at the Ruptura GitHub repository. It can be seen that the configuration is done by first creating a ``Components`` object where the user sets the molar fractions and the isotherms with key-value pairs. The ``isotherms`` field/key holds a list with two lists inside indicating a Dual-Site model. Each site, in this case, is modeled as a Langmuir-Freundlich isotherm.    

``` python
components = ruptura.Components([{
    "MoleculeName": "o-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 4.8, 1.28e-5, 1.4], ["Langmuir-Freundlich", 1.5, 1.6e-4, 0.7]]
}, {
    "MoleculeName": "m-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 4.8, 3.19e-5, 1.27], ["Langmuir-Freundlich", 1.2, 4.72e-5, 0.7]]
}, {
    "MoleculeName": "p-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 4.5, 2.3e-6, 1.7], ["Langmuir-Freundlich", 1.6, 1.46e-4, 0.7]]
}, {
    "MoleculeName": "ethylbenzene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 4.5, 7.38e-6, 1.5], ["Langmuir-Freundlich", 1.4, 1.6e-4, 0.7]]
}])

mix = ruptura.MixturePrediction(components=components,
                                DisplayName="CoBDP",
                                Temperature=443.0,
                                PressureStart=1e2,
                                PressureEnd=1e6,
                                NumberOfPressurePoints=1,
                                PressureScale="log")

data = mix.compute()
```

### PyIAST Syntax

### Julia Implementation

Running the implementations in Ruptura requires installation of visual studio build tools. To avoid problems in building the documentation of the Julia package, we will simply read the Ruptura's simulated results saved as a file the in this repository.

```@example ruptura 
using Langmuir
o_xylene     = MultiSite(LangmuirFreundlich(4.80, 1.28e-5, -0.0, 1.4,  -0.0), LangmuirFreundlich(1.5, 1.60e-4,  -0.0, 0.7, -0.0))
m_xylene     = MultiSite(LangmuirFreundlich(4.80, 3.19e-5, -0.0, 1.27, -0.0), LangmuirFreundlich(1.2, 4.72e-5,  -0.0, 0.7, -0.0))
p_xylene     = MultiSite(LangmuirFreundlich(4.50, 2.3e-6,  -0.0, 1.7,  -0.0), LangmuirFreundlich(1.6, 1.46e-4,  -0.0, 0.7, -0.0))
ethylbenzene = MultiSite(LangmuirFreundlich(4.50, 7.38e-6, -0.0, 1.5,  -0.0), LangmuirFreundlich(1.4, 1.60e-4,  -0.0, 0.7, -0.0))
models = [o_xylene, m_xylene, p_xylene, ethylbenzene]


ncomponents = 4
y = ones(ncomponents)/ncomponents
T = 443.0
p = 1e2
iast(models, p, T, y, FastIAS(), maxiters = 500, abstol = 1e-13, reltol = 1e-13)

```