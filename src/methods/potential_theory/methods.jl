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
    μ_ads = chemical_potential(eos, p_ads, T_ads, x_ads, :unknown)
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
        return (Δ, p_x) -> begin
            p = first(p_x)
            x = @view p_x[2:end]
            res_μ!(Δ, eos, p, T, x, μ, Ψ)
        end
    end
    return f
end

function solve_at_z(eos, p0_ads, T_ads, x0_ads, μ_bulk::M, Ψ::M, alg::A; verbose = true) where {M <: Number, A <: ChemPotentialMethod}
    abstol = alg.abstol
    reltol = alg.reltol
    f = create_res_func(eos, T_ads, μ_bulk, Ψ)
    P_result = Roots.find_zero(x -> to_newton(f, x), p0_ads, Roots.Newton(), abstol = abstol, reltol = reltol, verbose = verbose)
    _1 = one(T_ads)
    ρ = _1/volume(eos, P_result, T_ads, x0_ads, :stable)
    
    # For single component, always return true for convergence (Roots.jl will throw if it doesn't converge)
    return ρ, P_result, true
end

function solve_at_z(eos, p0_ads, T_ads, x0_ads::M, μ_bulk::M, Ψ::M, alg::A; verbose = true) where {M <: AbstractVector, A <: ChemPotentialMethod}
    f! = create_res_func(eos, T_ads, μ_bulk, Ψ)
    x0 = [p0_ads; x0_ads...]
    abstol = alg.abstol
    reltol = alg.reltol
    options = NEqOptions(f_abstol = abstol, f_reltol = reltol, maxiter = 10_000)
    sol = nlsolve(f!, x0, options = options)
    _1 = one(T_ads)
    p, x = sol.info.solution[1], sol.info.solution[2:end]
    ρ = _1/volume(eos, p, T_ads, x, :unknown).*x
    
    # Check convergence: norm of residual should be small
    res_norm = norm(sol.info.best_residual)
    converged = res_norm < abstol || res_norm < reltol * norm(sol.info.ρF0)
    
    # Return density, pressure, and convergence flag
    return ρ, p, converged
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
        ρᵢ, Pᵢ, converged = solve_at_z(eos, Pᵢ, T, xᵢ, μ_bulk, Ψ, alg, verbose = false)
        
        # Warn if not converged - indicates potential issue with solution at this point
        # Note: ChemPotentialMethod can converge to incorrect solutions at very high potentials
        # (near the wall) when using previous point as initial guess. FugacityCoefficientMethod
        # is more robust for such cases.
        if !converged && verbose
            @warn "solve_at_z did not converge at i=$i (z=$z, Ψ=$Ψ). Solution at this point may be inaccurate."
        end
        
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


#FUGACITY COEFFICIENT METHOD FOR POTENTIAL THEORY

mutable struct FugacityCoefficientMethod{T <: Real}
    abstol::T
    reltol::T
    maxiter_outer::Int
    maxiter_inner::Int
    x0::Union{PTASolution, Nothing}
    detect_phase_transition::Bool
end

function FugacityCoefficientMethod(; abstol=1e-7, reltol=1e-7, 
                                    maxiter_outer=150, maxiter_inner=100,
                                    x0=nothing, detect_phase_transition=false)
    FugacityCoefficientMethod(abstol, reltol, maxiter_outer, maxiter_inner, 
                              x0, detect_phase_transition)
end

function FugacityCoefficientMethod(prob::PTAProblem; kwargs...)
    x0 = PTA_x0(prob)
    FugacityCoefficientMethod(; x0=x0, kwargs...)
end

"""
    solve_at_potential_fugacity(eos, P_init, z_init, T, z_bulk, f_bulk, P_bulk, ε_vec, alg; verbose=false)

Solve for P(ε) and z(ε) at a single potential value using fugacity coefficient method
with stability checks to avoid oscillation.

Paper equations:
- Composition: z^i(ε) = [z_g^i * f_g^i * P_g * exp(ε^i/RT)] / [f^i(P,z)]
- Constraint:  Σ_i [z_g^i * f_g^i * P_g * exp(ε^i/RT)] / [f^i(P,z)] - 1 = 0

Uses Clapeyron's lnϕ and ForwardDiff for automatic differentiation.

Algorithm (with stability check from Clapeyron bubble pressure):
OUTER LOOP (Newton on pressure):
  1. Save checkpoint: z_restart, P_restart
  2. INNER LOOP (successive substitution on composition):
     a. Calculate f(P, z) using lnϕ
     b. Update z from fugacity equation
     c. Check convergence on z
     d. STABILITY CHECK: ||z - z_bulk|| < tol_stability?
        - If yes: restore checkpoint, exit inner loop, force P update
        - If no: continue
  3. Calculate pressure residual and derivative using ForwardDiff
  4. Newton update on P with damping
  5. Check convergence
"""
function solve_at_potential_fugacity(eos, P_init, z_init, T, z_bulk, f_bulk, P_bulk, ε_vec, alg; verbose=false)
    
    # Unpack
    RT = Rg(eos) * T
    n_comp = length(z_bulk)
    
    # Pre-compute constant term: C^i = z_g^i * f_g^i  * exp(ε^i/RT)
    C = z_bulk .* f_bulk .* exp.(ε_vec ./ RT)
    
    # Tolerances
    tol_z = alg.abstol
    tol_P = alg.reltol
    tol_stability = abs2(cbrt(tol_z))  # Stability threshold (Clapeyron pattern)
    
    # Initialize
    P = P_init
    z = copy(z_init)
    
    # OUTER LOOP: Newton iteration on pressure
    for j in 1:alg.maxiter_outer
        
        # Checkpoint for stability restart
        z_restart = copy(z)
        P_restart = P
        
        valid_iter = true  # Flag for valid inner loop convergence
        
        # INNER LOOP: Successive substitution on composition
        for i in 1:alg.maxiter_inner
            
            # Calculate fugacity coefficients at current (P, z)
            # Using Clapeyron's lnϕ function
            lnϕ_z = lnϕ(eos, P, T, z)[1]
            ϕ_z = exp.(lnϕ_z)
            f_z = ϕ_z .* P  # Fugacities
            
            # Update composition: z^i = C^i / f^i(P,z)
            z_new = C ./ f_z
            z_new = z_new / sum(z_new) # Normalize
            
            # Check convergence on composition
            if norm(z_new - z, Inf) < tol_z
                z = z_new
                break
            end
            
            #═════════════════════════════════════════════════════
            # STABILITY CHECK (prevents oscillation/trivial solution)
            #═════════════════════════════════════════════════════
            # K-value based stability: sum of squared log deviations from bulk
            # Non-allocating implementation from Clapeyron
            K_stability = zero(Base.promote_eltype(z_new, z_bulk))
            for k in eachindex(z_new)
                K_stability += log(z_new[k] / z_bulk[k])^2
            end
            
            if K_stability < tol_stability
                # Segregated composition too similar to bulk
                # → Iteration going wrong (weak field, wrong P, near critical)
                
                verbose && @warn "Stability check triggered at P=$P: K_stability=$K_stability < $tol_stability"
                
                # Restore to checkpoint
                z = z_restart
                P = P_restart
                
                # Mark iteration as invalid (skip derivative calc, force P update)
                valid_iter = false
                
                # Exit inner loop
                break
            end
            #═════════════════════════════════════════════════════
            
            # Update for next inner iteration
            z = z_new
        end
        
        # If stability check triggered, skip to pressure update
        if !valid_iter
            verbose && @info "Forcing pressure update due to stability check"
        end
        
        # Calculate fugacity coefficients and residual
        lnϕᵢ = lnϕ(eos, P, T, z)[1]
        ϕᵢ = exp.(lnϕᵢ)
        fᵢ = ϕᵢ .* P
        
        # Residual: F = Σ[Cⁱ/fⁱ] - 1
        F = sum(C ./ fᵢ) - 1.0
        
        # Analytical derivative: ∂F/∂P = -Σᵢ[Cⁱ/fⁱ²]·∂fⁱ/∂P
        # where ∂fⁱ/∂P = ϕⁱ(1 + P·∂ln(ϕⁱ)/∂P)
        dlnϕdP = ∂lnϕ∂n∂P(eos, P, T, z)[3]  # Optimized function
        dfᵢdP = ϕᵢ .* (1.0 .+ P .* dlnϕdP)
        ∂F∂P = -sum((C ./ (fᵢ .^ 2)) .* dfᵢdP)
        
        # Newton update with damping
        ΔP = F / ∂F∂P
        P_new = P - clamp(ΔP, -0.4*P, 0.4*P)  # Limit step size to ±40%
        
        # Check convergence
        if valid_iter && (abs(F) < tol_P && abs(ΔP/P) < tol_P)
            verbose && @info "Converged at P=$P_new, ||z-z_bulk||=$(norm(z-z_bulk, Inf))"
            return P_new, z, true
        end
        
        P = P_new
    end
    
    # Did not converge
    @warn "solve_at_potential_fugacity did not converge after $(alg.maxiter_outer) iterations"
    return P, z, false
end

"""
    solve_PTAProblem(prob::PTAProblem, alg::FugacityCoefficientMethod; verbose=true)

Main solver using fugacity coefficient method. Loops over potential values ε (not distance z).

At each potential value ε, solves the coupled system:
- Composition equation: z^i(ε) = [z_g^i * f_g^i * P_g * exp(ε^i/RT)] / [f^i(P(ε),z(ε))]
- Summation constraint: Σ_i z^i(ε) = 1 (implicitly via pressure equation)

Returns PTASolution with density and composition profiles.

Note: Only valid for multicomponent systems. For single components, use ChemPotentialMethod.
"""
function solve_PTAProblem(prob::PTAProblem, alg::FugacityCoefficientMethod; verbose=true)
    
    eos = prob.system.model
    potential_model = prob.system.structure.potential
    T = prob.T
    P_bulk = prob.P
    z_bulk = prob.x
    
    # Check for multicomponent
    if length(z_bulk) == 1
        error("FugacityCoefficientMethod only applies to multicomponent systems. Use ChemPotentialMethod for single components.")
    end
    
    # Compute bulk fugacity coefficients using Clapeyron
    lnϕ_bulk = lnϕ(eos, P_bulk, T, z_bulk)[1]
    ϕ_bulk = exp.(lnϕ_bulk)
    f_bulk = ϕ_bulk .* P_bulk  # Bulk fugacities
    
    # Initialize solution (already contains grid)
    sol = isnothing(alg.x0) ? PTA_x0(prob) : alg.x0


    # Initial guess: start from bulk conditions
    P_ε = P_bulk
    z_ε = copy(z_bulk)
    
    # Loop over potential values (from high to low, i.e., z0 to boundary)
    for i in reverse(eachindex(sol.z))
        ε = potential(potential_model, sol.z[i])
        
        # For MultiComponentDRA, ε is a vector (one per component)
        # For DRA (single potential), ε is a scalar - broadcast to all components
        ε_vec = isa(ε, AbstractVector) ? ε : fill(ε, length(z_bulk))
        
        # Solve at this potential value
        P_ε, z_ε, converged = solve_at_potential_fugacity(
            eos, P_ε, z_ε, T, z_bulk, f_bulk, P_bulk, ε_vec, alg, 
            verbose=verbose
        )
        
        if !converged
            @warn "Failed to converge at ε=$ε (i=$i)"
        end
        
        # Calculate density and store results
        v_ε = volume(eos, P_ε, T, z_ε, :stable)
        ρ_total = 1.0 / v_ε  # Total molar density
        ρ_comp = ρ_total .* z_ε  # Component densities
        
        sol.P[i, 1] = P_ε
        sol.x[i, :] .= z_ε
        sol.ρ[i, :] .= ρ_comp
    end
    
    sol.retcode = :success
    alg.x0 = sol
    
    return sol
end

function loading(prob::PTAProblem, alg::FugacityCoefficientMethod; verbose=false)
    sol = solve_PTAProblem(prob, alg; verbose=verbose)
    x_bulk = prob.x
    ρ_bulk = prob.bulkcondition[1]
    return loading(sol, ρ_bulk, x_bulk)
end

export PTAProblem, PTASolution, ChemPotentialMethod, FugacityCoefficientMethod
