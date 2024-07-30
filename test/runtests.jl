using AdsorbedSolutionTheory
using Test

#we test that definitions of loading and sp_res are consistent.
function test_sp_res_loading(model,prange,T)
    for p in prange
        n = loading(model,p,T)
        n₀ = loading_ad(model,p,T)
        Π = sp_res(model,p,T)
        Π₀ = sp_res_numerical(model,Π,T)
        @test n ≈ n₀
        @test Π ≈ Π₀
    end
end

@testset "AdsorbedSolutionTheory.jl" begin
    # Write your tests here.
end
