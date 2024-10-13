# Models in adsorption equilibrium

Adsorption is a surface phenomenon where atoms, ions, or molecules from a gas, liquid, or dissolved solid adhere to the surface of another material. This process leads to the formation of a thin film of the adsorbed substance, known as the adsorbate, on the surface of the material it adheres to, called the adsorbent.

In thermodynamics, adsorption models play a similar role to equations of state in fluid systems, describing the equilibrium properties of adsorption processes. These models are crucial for understanding and predicting the thermodynamics of adsorption in both single and multicomponent systems.

## Single-component adsorption

Over the last decade, three fundamental approaches have led to the development of a wide range of adsorption isotherm models, including well-known ones such as Langmuir, Freundlich, Dubinin-Radushkevich, Temkin, Toth, and many others.

The first approach focuses on kinetics, where adsorption and desorption rates are equal, establishing adsorption equilibrium as a dynamic process. The second approach is rooted in thermodynamics, offering a basis for deriving various forms of adsorption isotherm models. The third approach emphasizes the generation of characteristic curves to describe adsorption behavior.

### Isosteric heat of adsorption

The heat of adsorption is a critical design parameter in adsorptive gas separation units. During adsorption, heat is released as adsorbate molecules transition to a lower energy state on the surface of the adsorbent compared to their higher energy state in the bulk gas phase. This exothermic process significantly impacts both the efficiency and operational conditions of adsorption systems. For a single component, the isosteric heat is given by:

$Q_{st} = -T*(V_g - V_a)*\left( \frac{dP_i}{dT} \right)\rvert_{(N_i,A)}$ 

where $Q_st$ is the isosteric heat of the component being adsorbed, $T$ is the temperature, $V_g$ is the molar volume of the component in gas phase, $V_a$ is the molar volume of the component in adsorbed phase, $N_i$ is the amount of component adsorbed of the component

When the isotherm is of the form $N_i = f(T, P_i)$, one can write:

$Q_{st, i} = -T*(V_g - V_a)*\left( \frac{\frac{\partial N_i}{\partial T}\rvert_P}{\frac{\partial N_i}{\partial P}\rvert_T} \right)$ 


## Multi component adsorption

The basic equations of the IAST are the analogue of Raoult's law in vapour–liquid equilibrium:

$Py_i = P_i^0(\pi)x_i$ 

where

$\pi = \pi_i = \int_{0}^{P_i^0} \frac{N_i^0(P)}{P}dP$ for $i = 1,...,N_c$ 

$\sum_i^{N_c} x_i = 1$ 


Combining (1) and (3), the following nonlinear solve is set to:

$f(\pi) = 1 - \sum_1^{N_c}\frac{Py_i}{P_i^0\left(\pi\right)}$ = 0 













