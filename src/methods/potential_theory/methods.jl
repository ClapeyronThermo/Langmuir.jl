struct PTAProblem{S <: PTASystem, TP <: Real, CondType, XType <: AbstractVector{<: TP}}
    system::S
    T::TP
    P::TP
    bulkcondition::CondType
    x::XType

    function PTAProblem(system::S, T, P, bulkcondition, x) where {S <: PTASystem}
        x = size_check(system.model, x)
        T, P = promote(T, P)
        typeX = typeof(x)
        typeT = typeof(T)
        typebulk = typeof(bulkcondition)
        new{S, typeT, typebulk, typeX}(system, T, P, bulkcondition, x)
    end
end

function size_check(model, x)
    
    if length(model) == 1 && isa(x, Real)
        x = [x]
    end

    length(x) == length(model) || throw(ArgumentError("x length must match model length"))

    return x
end

function PTAProblem(T, P, x; eos, potential::V) where {V <: AbstractPotential}
    system = PTASystem(eos, potential)
    v = volume(eos, P, T, x, :stable)
    #is_stable = VT_isstable(eos, v, T, x)
    #is_stable || error("Phase is not stable")
    bulk_μ = ifelse(length(x) == 1, first(VT_chemical_potential(eos, v, T, x)), VT_chemical_potential(eos, v, T, x))
    n = sum(x)
    return PTAProblem(system, T, P, (n/v, bulk_μ), x)
end

mutable struct PTASolution{PType,XType,RType,ZType, S <: Symbol}
    P::PType
    x::XType
    ρ::RType
    z::ZType
    retcode::S

    function PTASolution(P, x, ρ, z, retcode = :unconverged)
        _P = typeof(P)
        _X = typeof(x)
        _R = typeof(ρ)
        _Z = typeof(z)
        _S = typeof(retcode)
        return new{_P, _X, _R, _Z, _S}(P, x, ρ, z, retcode)
    end
end

#= mutable struct PTAIteration{PType, XType, RzType}
    P::PType
    x::XType
    ρ::RzType
    function PTAIteration(P, x, ρ)
        Pz = typeof(P)
        Xz = typeof(x)
        Rz = typeof(ρ)
        return new{Pz, Xz, Rz}(P, x, ρ)
    end
end =#

function initialize_solution(prob::PTAProblem)
    grid_size = prob.system.structure.ngrids[1]
    ncomps = length(prob.system.model)
    P = Array{typeof(prob.P)}(undef, grid_size, 1)
    x = Array{typeof(prob.P)}(undef, grid_size, ncomps)
    ρ = Array{typeof(prob.P)}(undef, grid_size, ncomps) # Initialize as a vector of vectors
    Pz = zero(typeof(prob.P))
    Xz = Vector{typeof(prob.P)}()
    ρz = Vector{typeof(prob.P)}()
    z = create_grid(prob.system.structure)
    return PTASolution(P, x, ρ, z, :unconverged)
end


mutable struct ChemPotentialMethod{T <: Real}
    abstol::T
    reltol::T
    x0::Union{PTASolution, Nothing}
end

function ChemPotentialMethod(; abstol=1e-7, reltol=1e-7, x0=nothing)
    ChemPotentialMethod(abstol, reltol, x0)
end

# Constructor that takes a PTAProblem and sets x0 using PTA_x0(prob)
function ChemPotentialMethod(prob::PTAProblem; abstol=1e-7, reltol=1e-7)
    x0 = PTA_x0(prob)
    ChemPotentialMethod(abstol, reltol, x0)
end

function PTA_x0(prob::PTAProblem)
    sol = initialize_solution(prob)
    eosmodel = prob.system.model
    potential_model = prob.system.structure.potential
    P = prob.P
    T = prob.T
    x = prob.x
    ρ = molar_density(idealmodel(eosmodel), P, T, x)
    RgT⁻¹ = one(eltype(P))/(Rg(eosmodel) * T)

    for (i, z) ∈ enumerate(sol.z)
        ρᵢ = ρ .* x .* exp.(potential(potential_model, z).*RgT⁻¹)
        ∑ρ = sum(ρᵢ)
        P = Rg(eosmodel) * T * ∑ρ
        x = ρᵢ ./ ∑ρ
        sol.P[i, 1] = P
        sol.x[i, :] .= x
        sol.ρ[i, :] .= ρᵢ
    end

    return sol
end


function res_μ(eos, p_ads, T_ads, x_ads::M, μ_bulk::M, Ψ::M) where {M <: Number}
    v_ads = volume(eos, p_ads, T_ads, x_ads)
    #is_stable = VT_isstable(eos, v_ads, T_ads, x_ads)
    #is_stable || @warn "Adsorbed phase is not stable - results may be inaccurate"
    μ_ads = VT_chemical_potential(eos, v_ads, T_ads, x_ads)
    Δ = (μ_bulk .+ Ψ) .- μ_ads  
    return Δ
end

function res_μ!(Δ, eos, p_ads, T_ads, x_ads, μ_bulk, Ψ)
    #is_stable = isstable(eos, p_ads, T_ads, x_ads)
    #is_stable || @warn "Adsorbed phase is not stable - results may be inaccurate"
    μ_ads = chemical_potential(eos, p_ads, T_ads, x_ads, :stable)
    ∑x_ads = sum(x_ads)
    Δ .= [(μ_bulk .+ Ψ) .- μ_ads; ∑x_ads - 1.0]
    return Δ
end

# When only one component is present
function create_res_func(eos, T_ads, μ_bulk::M, Ψ::M) where {M <: Number}
    f = let eos = eos, T = T_ads, μ = μ_bulk, Ψ = Ψ
        return p -> res_μ(eos, p, T, 1.0, μ, Ψ)
    end
    return f
end

function create_res_func(eos, T_ads, μ_bulk::M, Ψ::M) where {M <: AbstractVector}
    #Δ = similar(μ_bulk, length(μ_bulk) + 1)
    f = let eos = eos, T = T_ads, μ = μ_bulk, Ψ = Ψ
        return (Δ, p_x) -> res_μ!(Δ, eos, first(p_x), T, p_x[2:end], μ, Ψ)
    end
    return f
end

function solve_at_z(eos, p0_ads, T_ads, x0_ads, μ_bulk::M, Ψ::M, alg::A; verbose = true) where {M <: Number, A <: ChemPotentialMethod}
    abstol = alg.abstol
    reltol = alg.reltol
    f = create_res_func(eos, T_ads, μ_bulk, Ψ)
    P_result = Roots.find_zero(x -> to_newton(f, x), p0_ads, Roots.Newton(), abstol = abstol, reltol = reltol, verbose = verbose)
    _1 = one(T_ads)
    return _1/volume(eos, P_result, T_ads, x0_ads, :stable), P_result
end

function solve_at_z(eos, p0_ads, T_ads, x0_ads::M, μ_bulk::M, Ψ::M, alg::A; verbose = true) where {M <: AbstractVector, A <: ChemPotentialMethod}
    f! = create_res_func(eos, T_ads, μ_bulk, Ψ)
    x0 = [p0_ads; x0_ads...]
    δ = similar(x0)
    f!(δ, x0)
    abstol = alg.abstol
    reltol = alg.reltol
    options = NEqOptions(f_abstol = abstol, f_reltol = reltol)
    sol = nlsolve(f!, x0, options = options)
    _1 = one(T_ads)
    p, x = sol.info.solution[1], sol.info.solution[2:end]
    return _1/volume(eos, p, T_ads, x, :stable).*x, p
end

function solve_PTAProblem(prob::PR, alg::A; verbose = true) where {PR <: PTAProblem, A <: ChemPotentialMethod}
    eos = prob.system.model
    _potential = prob.system.structure.potential
    T = prob.T
    P = prob.P
    x = prob.x
    μ_bulk = prob.bulkcondition[2]
    sol = isnothing(alg.x0) ? PTA_x0(prob) : alg.x0
    Pᵢ = P
    xᵢ = x

    for i ∈ reverse(eachindex(sol.z))
        z = sol.z[i]
        Ψ = potential(_potential, z)
        ρᵢ, Pᵢ = solve_at_z(eos, Pᵢ, T, xᵢ, μ_bulk, Ψ, alg, verbose = verbose)
        xᵢ = ρᵢ ./ sum(ρᵢ)
        sol.ρ[i, :] .= ρᵢ
        sol.P[i, 1] = Pᵢ
        sol.x[i, :] .= xᵢ
    end

    sol.retcode = :success
    alg.x0 = sol # update the algorithm's initial guess for future calls

    return sol
end

function loading(prob; solver = ChemPotentialMethod(prob))
    sol = solve_PTAProblem(prob, solver; verbose = false)
    x_bulk = prob.x
    return loading(sol, prob.bulkcondition[1], x_bulk)
end

function loading(sol::S, ρ_bulk, x_bulk; integrator = SimpsonsRule()) where {S <: PTASolution}
    yᵢ = sol.ρ .- ρ_bulk*x_bulk'
    z = sol.z
    problem = SampledIntegralProblem(yᵢ, z, dim = 1)
    solution = Integrals.solve(problem, integrator)
    return solution.u
end

export PTAProblem, PTASolution, ChemPotentialMethod
