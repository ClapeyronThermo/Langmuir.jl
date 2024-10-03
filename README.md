[![Build Status](https://github.com/ClapeyronThermo/Langmuir.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ClapeyronThermo/Langmuir.jl/actions/workflows/CI.yml?query=branch%3Amain) [![codecov](https://codecov.io/gh/ClapeyronThermo/Langmuir.jl/branch/main/graph/badge.svg?token=ZVGGR4AAFB)](https://codecov.io/gh/ClapeyronThermo/Langmuir.jl) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://clapeyronthermo.github.io/Langmuir.jl/dev)

![logo](/docs/Langmuir_logo.svg)

This package implements single adsorption and multicomponent adsorption through Ideal Absorbed Solution Theory (IAST). Some major features are:

    - Explicit single component adsorption properties: loading, bulk phase properties
    - Multicomponent adsorption 

## Examples:

```julia
using Langmuir

#example from doi.org/10.1002/aic.14684

v = @MultiSite{LangmuirS1,LangmuirS1} #we create a multisite model, consisting

    x1 = [1.468
        0.024
        0
        7.891
        0.001645
        0]

    x2 = [2.847
        0.028
        0.
        2.223
        1.228
        0.]

    x3 = [2.581
        0.84
        0.0
        2.901
        0.021
        0.0]

    #creation of isotherms from vectors or other iterables.
    #you can also use Langmuir.from_vec! or Langmuir.from_vec
    #to create (or fill) vectors from isotherm models
    m1,m2,m3 = Langmuir.from_vec(v,x1),Langmuir.from_vec(v,x2),AST.from_vec(v,x3)
    models = (m1,m2,m3)

    #calculate loading of a single isotherm
    #there are also functions to calculate the reduced spreading pressure, and inverse algorithms
    l1 = loading(m1,p,T)
    kh1 = henry_coefficient(m1,T)
    lmax = saturated_loading(m1,T)

    y = [0.5,0.25,0.25]
    T = 300
    p = 1000

    
    #total adsorbed amount, fractions of adsorbed components.
    q_tot,x,status = iast(models,p,T,y) 
```
