# 1. Traction-Free Elliptical Hole.jl
include("Common.jl")
using .Common

# ==========================================
# Analytical Solution (Inglis / Muskhelishvili)
# ==========================================

"""
    φ(ζ) for Traction Free Hole
    
    Γ = -(P1+P2)/4
    Γ' = (P1-P2)/2 * exp(-2iΛ)
    
    Formula: φ(ζ) = ΓRζ - (mΓ + conj(Γ')) * (R/ζ)
"""
function φ_traction_free(ζ, cav, mat, stress)
    P1, P2, Λ = stress.P1, stress.P2, stress.Λ
    R, m = cav.R, cav.m
    
    Γ  = -(P1 + P2) / 4.0
    Γp = (P1 - P2) / 2.0 * exp(-2im * Λ)
    
    term1 = Γ * R * ζ
    term2 = - (m * Γ + conj(Γp)) * (R / ζ)
    
    return term1 + term2
end

"""
    ψ(ζ) for Traction Free Hole
    
    Derived from boundary condition: φ + ω/ω' * conj(φ') + conj(ψ) = 0
    
    Formula: ψ(ζ) = Γ'Rζ - (conj(Γ) + m*conj(Γ')) * (R/ζ) - Term_Interaction
    Term_Interaction = (ζ * (1 + m*ζ^2) / (ζ^2 - m)) * φ'(ζ)
"""
function ψ_traction_free(ζ, cav, mat, stress)
    P1, P2, Λ = stress.P1, stress.P2, stress.Λ
    R, m = cav.R, cav.m
    
    Γ  = -(P1 + P2) / 4.0
    Γp = (P1 - P2) / 2.0 * exp(-2im * Λ)
    
    # Base terms
    term1 = conj(Γp) * R * ζ # Note: ψ usually has Γ' at infinity, but in conj form logic it varies.
    # Standard Muskhelishvili form for ψ usually has Γ'Rζ at infinity.
    # Let's use the explicit derived form for the hole:
    
    term1 = Γp * R * ζ
    term2 = - (conj(Γ) + m * conj(Γp)) * (R / ζ)
    
    # Interaction term: - [ω_bar(1/ζ) / ω'(ζ)] * φ'(ζ)
    # The term [ω_bar(1/ζ) / ω'(ζ)] simplifies to ζ(1+mζ²)/(ζ²-m)
    
    # We need φ'(ζ) analytically here for precision
    dφ = Γ * R + (m * Γ + conj(Γp)) * (R / (ζ^2))
    
    interaction = - (ζ * (1 + m * ζ^2) / (ζ^2 - m)) * dφ
    
    return term1 + term2 + interaction
end

# ==========================================
# Simulation Setup
# ==========================================

function run_simulation()
    # 1. Geometry: Ellipse 2:1 ratio
    a0 = 2.0
    b0 = 1.0
    cav = Common.EllipticalCavity(a0, b0)
    
    # 2. Material (Steel-ish)
    mat = Common.ElasticMaterial(E=200e9, ν=0.3)
    
    # 3. Load: Uniaxial Tension along Y-axis (P2 > 0, P1 = 0)
    # This creates high stress concentration at x=a (the tips)
    P_inf = 100.0 # MPa
    stress = Common.FarFieldStress(0.0, P_inf, 0.0) 
    
    println("Simulating Traction-Free Elliptical Hole...")
    println("Load: Uniaxial Tension $(P_inf) MPa along Y-axis")
    println("Geometry: a=$(a0), b=$(b0)")
    
    # 4. Plot
    Common.plot_results(cav, mat, stress, φ_traction_free, ψ_traction_free)
end

run_simulation()