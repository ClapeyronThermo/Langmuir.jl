using Langmuir
import Langmuir: loading_ad, sp_res, to_vec, sp_res_numerical, isosteric_heat, Rgas, from_vec, fit, pressure, temperature, x0_guess_fit
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
        ΔH = isosteric_heat(cL_test, 101325.0, 298.) ≈ -cL_test.E
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
end
