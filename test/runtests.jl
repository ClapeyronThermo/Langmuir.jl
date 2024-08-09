using AdsorbedSolutionTheory
import AdsorbedSolutionTheory: loading_ad, sp_res, sp_res_numerical, isosteric_heat, Rgas
using Test

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
        cL_test = Langmuir(1.727, 16.71e-10, -16152.50)

        # Loading
        test_sp_res_loading(cL_test, [1e5, 2e5, 3e5], 298.)

        # Isosteric heat
        V = Rgas(cL_test)*298.0/101325.0
        ΔH = isosteric_heat(cL_test, V, 101325.0, 298.) ≈ -cL_test.E
    end
end

@testset "IAST" begin
    #10.1002/aic.14684, S.I
    v = @MultiSite{Langmuir,Langmuir}
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

end