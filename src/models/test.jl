using Langmuir, Plots

P = 1000
T = 283

modelL = LangmuirS1(6.84,1.15e-12,-54000.)
LL = loading(modelL,P,T)

modelI1 = Ising(12.3,2.35e-12,-43900.,1.04e-10,-20500.)
LI1 = loading(modelI1,P,T)

modelI2 = Ising(10.4,5.17e-12,-46800.,2.34e-7,-1330.)
LI2 = loading(modelI2,P,T)


combo = MultiSite(modelL,modelI1,modelI2)
LTOT = loading(combo1,P,T)


pressures = 0:250:4000

loadings0 = [loading(combo,pressure,243) for pressure in pressures]
loadings1 = [loading(combo,pressure,263) for pressure in pressures]
loadings2 = [loading(combo,pressure,283) for pressure in pressures]
loadings3 = [loading(combo,pressure,303) for pressure in pressures]
loadings4 = [loading(combo,pressure,323) for pressure in pressures]
loadings5 = [loading(combo,pressure,343) for pressure in pressures]

plot(pressures, loadings1,xlabel="Pressure(Pa)",ylabel="loading(mol⋅kg⁻¹)", label="263K", lw=2,markershape=:circle)
plot!(pressures, loadings0, label="243K", lw=2,markershape=:circle)
plot!(pressures, loadings2, label="328K", lw=2,markershape=:circle)
plot!(pressures, loadings3, label="303K", lw=2,markershape=:circle)
plot!(pressures, loadings4, label="323K", lw=2,markershape=:circle)
plot!(pressures, loadings5, label="343K", lw=2,markershape=:circle)

