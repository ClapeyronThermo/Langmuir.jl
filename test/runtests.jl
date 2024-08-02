using AdsorbedSolutionTheory
import AdsorbedSolutionTheory: loading_ad, sp_res, sp_res_numerical
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

@testset "AdsorbedSolutionTheory.jl" begin
    #Langmuir testing
    cL_test = Langmuir(1.727, 16.71e-10, -16152.50)
    test_sp_res_loading(cL_test, [1e5, 2e5, 3e5], 298.)
end


