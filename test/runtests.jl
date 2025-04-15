using Langmuir
import Langmuir: loading_ad, sp_res, to_vec, sp_res_numerical, isosteric_heat, Rgas, from_vec, fit, pressure, temperature, x0_guess_fit
import Langmuir: gibbs_excess_free_energy, activity_coefficient
import Langmuir: IsothermFittingProblem, DEIsothermFittingSolver
using Test
const LG = Langmuir



#we test that definitions of loading and sp_res are consistent.

function test_sp_res_loading(model, prange, T)

    for p in prange
        n = loading(model,p,T)
        n₀ = loading_ad(model,p,T)
        Π = sp_res(model,p,T)
        Π₀ = sp_res_numerical(model, p,T)
        @test n ≈ n₀
        @test Π ≈ Π₀
    end
end

@testset "isotherm consistency" begin
    #Langmuir testing
    @testset "langmuir" begin
        cL_test = LangmuirS1(1.727, 16.71e-10, -16152.50)

        # Loading
        test_sp_res_loading(cL_test, [1e5, 2e5, 3e5], 298.)

        # Isosteric heat
        ΔH = isosteric_heat(cL_test, 101325.0, 298.) ≈ cL_test.E
    end
end

@testset "isotherm fitting" begin
    @testset "langmuir" begin
        p = range(0.0, 101325.0*1.0, length = 50) |> collect
        t = range(31.6 + 273.15, 50.0 + 273.15, length = 5) |> collect
        P = vec(p'.*ones(length(t)))
        T = vec(ones(length(p))' .* t)
        pt_ = hcat(P, T)
        lang = LangmuirS1(5.073, 0.67e-10, -39300.55246960819)
        σ = 0.05
        l = map(pT -> abs(loading(lang, pT[1], pT[2]) + σ*randn()), eachrow(pt_))
        table = (;P,l,T)
        d = isotherm_data(table, :P, :l, :T)

        prob = IsothermFittingProblem(LangmuirS1, d, abs2)
        alg = DEIsothermFittingSolver(logspace = true)
        loss_fit, fitted_isotherm = fit(prob, alg)

        @test (abs(sqrt(loss_fit) - σ)/σ)*100.0 < 10.0 #relative error smaller than 5% 
    end

    @testset "quadratic" begin
        p = range(0.0, 101325.0*1.0, length = 50) |> collect
        t = range(31.6 + 273.15, 50.0 + 273.15, length = 5) |> collect
        P = vec(p'.*ones(length(t)))
        T = vec(ones(length(p))' .* t)
        pt_ = hcat(P, T)
        quad_ = Quadratic(0.67e-10, 0.37e-11, 4.073, -37300.9, -23300.55)
        σ = 0.05
        l = map(pT -> abs(loading(quad_, pT[1], pT[2]) + σ*randn()), eachrow(pt_))
        table = (;P,l,T)
        d = isotherm_data(table, :P, :l, :T)

        x0 = to_vec(x0_guess_fit(Quadratic, d))
        lb = (1e-35, 1e-35, 1e-29, -5_000., -5_000.)
        ub = (1e-3, 1e-3, 100., -80_000., -80_000.)
        prob = IsothermFittingProblem(Quadratic{eltype(d)}, d, nothing, abs2, x0, lb, ub) #Bounds have to be manually tweaked. Default interval is too large
        alg = DEIsothermFittingSolver(max_steps = 5000, logspace = true)
        loss_fit, fitted_isotherm = fit(prob, alg)

        @test (abs(sqrt(loss_fit) - σ)/σ)*100.0 < 10.0 #relative error smaller than 5% 
    end
end

@testset "ThermodynamicLangmuir" begin
    tlang = ThermodynamicLangmuir(2.468, 7.03e-10, -2540*Rgas(LangmuirS1{Float64}), -611.63)
    analyical_gammas(T, θ) = (exp(θ[2]^2*tlang.Bᵢᵩ/T*(exp(-0.3*tlang.Bᵢᵩ/T) - 1)/(θ[1]*exp(-0.3*tlang.Bᵢᵩ/T) + θ[2])^2),
    exp(θ[1]^2*-tlang.Bᵢᵩ/T*(exp(-0.3*-tlang.Bᵢᵩ/T) - 1)/(θ[1] + θ[2]*exp(-0.3*-tlang.Bᵢᵩ/T))^2)                                                                                         )
    @test analyical_gammas(273.15, [0.4, 0.6])[1] ≈ activity_coefficient(tlang, 273.15, [0.4, 0.6])[1]
    @test analyical_gammas(273.15, [0.4, 0.6])[2] ≈ activity_coefficient(tlang, 273.15, [0.4, 0.6])[2]
    
    tlang2 = ThermodynamicLangmuir(2.468, 7.03e-10, -2540*Rgas(LangmuirS1{Float64}), 0.0)
    @test loading(tlang2, 101325.0, 298.15) ≈ loading(LangmuirS1(tlang2.M, tlang2.K₀, tlang2.E), 101325.0, 298.15)
end

@testset "gibss free energy multicomponent" begin
    tl1 = ThermodynamicLangmuir(2.00, 7.0e-5, -10_000.0, -600.0)
    tl2 = ThermodynamicLangmuir(3.00, 7.0e-6, -20_000.0, -100.0)
    models = ThermodynamicLangmuirModels(tl1, tl2)
    nrtl = aNRTLModel(models)
    x = [0.8, 0.2]
    T = 300.0
    τ₁ = tl1.Bᵢᵩ/T
    τ₂ = tl2.Bᵢᵩ/T
    τ₁₂ = τ₁ - τ₂
    G₁₂ = exp(-0.3*(τ₁ - τ₂))
    analytical_Gᴱ = x[1]*x[2]*τ₁₂*(G₁₂ - 1.0)/(x[1]*G₁₂ + x[2])
    @test gibbs_excess_free_energy(nrtl, T, x) ≈ analytical_Gᴱ
end


@testset "Multicomponent Extended Langmuir" begin
    Lang1 = LangmuirS1(1.727, 16.71e-10, -16152.50)
    Lang2 = LangmuirS1(4.0, 6.71e-8, -14152.50)
    multilang = ExtendedLangmuir(Lang1, Lang2)
    y = [0.3, 0.7]
    p = 101325.0
    T = 400.0
    @test loading(multilang, p, T, y)[1] ≈ Lang1.M * Lang1.K₀*exp(-Lang1.E/(8.31446261815324*T)) * p * y[1] / (1.0 + 
    Lang2.K₀*exp(-Lang2.E/(8.31446261815324*T))*p*y[2] + Lang1.K₀*exp(-Lang1.E/(8.31446261815324*T))*p*y[1])
    @test loading(multilang, p, T, y)[2] ≈ Lang2.M * Lang2.K₀*exp(-Lang2.E/(8.31446261815324*T)) * p * y[2] / (1.0 +
    Lang2.K₀*exp(-Lang2.E/(8.31446261815324*T))*p*y[2] + Lang1.K₀*exp(-Lang1.E/(8.31446261815324*T))*p*y[1])
end


@testset "Multicomponent Extended Multisite Langmuir" begin
    Lang1 = MultiSite(LangmuirS1(1.727, 16.71e-10, -16152.50), LangmuirS1(4.0, 6.71e-8, -14152.50))
    Lang2 = MultiSite(LangmuirS1(2.727, 16.71e-12, 1.3*-16152.50), LangmuirS1(7.0, 6.71e-5, 2*-14152.50))
    multilang = ExtendedLangmuir(Lang1, Lang2)
    y = [0.6, 0.4]
    p = 2*101325.0
    T = 400.0

    den_highE = 1.0 + Lang2.isotherms[2].K₀*exp(-Lang2.isotherms[2].E/(8.31446261815324*T))*p*y[2] + 
              Lang1.isotherms[1].K₀*exp(-Lang1.isotherms[1].E/(8.31446261815324*T))*p*y[1]

    den_lowE = 1.0 + Lang2.isotherms[1].K₀*exp(-Lang2.isotherms[1].E/(8.31446261815324*T))*p*y[2] +
                Lang1.isotherms[2].K₀*exp(-Lang1.isotherms[2].E/(8.31446261815324*T))*p*y[1]

    @test loading(multilang, p, T, y)[1] ≈ Lang1.isotherms[1].M * Lang1.isotherms[1].K₀*exp(-Lang1.isotherms[1].E/(8.31446261815324*T)) * p * y[1] / den_highE + 
    Lang1.isotherms[2].M * Lang1.isotherms[2].K₀*exp(-Lang1.isotherms[2].E/(8.31446261815324*T)) * p * y[1] / den_lowE

    @test loading(multilang, p, T, y)[2] ≈ Lang2.isotherms[1].M * Lang2.isotherms[1].K₀*exp(-Lang2.isotherms[1].E/(8.31446261815324*T)) * p * y[2] / den_lowE + 
    Lang2.isotherms[2].M * Lang2.isotherms[2].K₀*exp(-Lang2.isotherms[2].E/(8.31446261815324*T)) * p * y[2] / den_highE
            
end


@testset "IAST" begin
    #10.1002/aic.14684, S.I
    v = @MultiSite{LangmuirS1,LangmuirS1}
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

    m1,m2,m3 = from_vec(v,x1),from_vec(v,x2),from_vec(v,x3)
    models = (m1,m2,m3)
    y = [0.5,0.25,0.25]
    T = 300
    p = 1000

    q_tot = 5.331793232381403
    x_res = [0.12076215703820294, 0.4075271543741945, 0.4717106885876027]
    x0 = [0.12,0.41,0.47]

    @test iast(models,p,T,y)[1] ≈ q_tot
    @test iast(models,p,T,y,x0 = x0)[1] ≈ q_tot
    @test iast(models,p,T,y,FastIAS())[1] ≈ q_tot
    @test iast(models,p,T,y,IASTNestedLoop())[1] ≈ q_tot
    @test iast(models,p,T,y,IASTNestedLoop(),x0 = x0)[1] ≈ q_tot

    ethane_isotherm = Quadratic{Float64}(2.5087482698420104e-7, 2.1373377526197646e-19, 3.641079631515442, -6898.20708031339, -47789.60001500269)
    ethylene_isotherm = Quadratic{Float64}(2.599227350906123e-8, 7.128313806215397e-19, 3.832139235999132, -11790.383728687304, -41702.74723166111)
    p2 = 101325.0
    T2 = 303.0
    models2 = (ethane_isotherm,ethylene_isotherm)
    yx = range(0.0, 1.00, 51)
    for yi in yx
        y = [yi,1-yi]
        n1,w1,status1 = iast(models2,p2,T2,y,FastIAS())
        n2,w2,status2 = iast(models2,p2,T2,y,IASTNestedLoop())
        @test status1 == :success
        @test status2 == :success
        @test n1 ≈ n2
        @test w1 ≈ w2
    end
end



