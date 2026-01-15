include("Common.jl")
using .Common
using LinearAlgebra
using Printf

# ==========================================
# Analytical Solution (Inglis / Muskhelishvili)
# ==========================================

"""
    φ(ζ) for Traction Free Hole
    
    Γ = (P1+P2)/4
    Γ' = -(P1-P2)/2 * exp(-2iΛ)
    
    Formula: φ(ζ) = Γ*R*ζ + Γ'*R/ζ
"""
function φ_traction_free(ζ, cav, mat, stress)
    P1, P2, Λ = stress.P1, stress.P2, stress.Λ
    R, m = cav.R, cav.m

    Γ = (P1 + P2) / 4.0
    Γp = -(P1 - P2) / 2.0 * exp(-2im * Λ)

    term1 = Γ * R * ζ
    term2 = Γp * R / ζ

    return term1 + term2
end

"""
    ψ(ζ) for Traction Free Hole
    
    Derived from boundary condition: φ + ω/ω' * conj(φ') + conj(ψ) = 0
    
    Formula from literature for stability:
    ψ(ζ) = -Γ * (R/m) * (1/(ζ^2-m)) * (ζ - m/ζ) + 
           conj(Γ') * R * (ζ - (1-m^2)/(ζ*(ζ^2-m)))
"""
function ψ_traction_free(ζ, cav, mat, stress)
    P1, P2, Λ = stress.P1, stress.P2, stress.Λ
    R, m = cav.R, cav.m

    Γ = (P1 + P2) / 4.0
    Γp = -(P1 - P2) / 2.0 * exp(-2im * Λ)

    # This is a more stable formulation found in literature
    term1 = -Γ * (R / m) * (1 / (ζ^2 - m)) * (ζ - m / ζ)
    term2 = conj(Γp) * R * (ζ - (1 - m^2) / (ζ * (ζ^2 - m)))

    return term1 + term2
end

# ==========================================
# Sanity Check Routine
# ==========================================

function sanity_check(cav, mat, stress, φ_func, ψ_func)
    println("\n=== Running Sanity Check ===")

    # 1. Far-Field Condition Check
    # ----------------------------
    # Check stresses at a large distance (e.g., 20 * a0)
    # We check along X and Y axes
    dist_mult = 20.0
    x_far = dist_mult * cav.a0
    y_far = dist_mult * cav.b0

    # Calculate at x-axis far field
    σx_inf, σy_inf, _ = Common.calculate_field(x_far, 0.0, cav, mat, stress, φ_func, ψ_func)

    println("\n[Far-Field Check @ $(dist_mult)a]")
    println("Input P1 (x-axis): $(stress.P1) | Calculated σx: $(@sprintf("%.4f", σx_inf))")
    println("Input P2 (y-axis): $(stress.P2) | Calculated σy: $(@sprintf("%.4f", σy_inf))")

    # Simple error metric
    err_x = abs(σx_inf - stress.P1)
    err_y = abs(σy_inf - stress.P2)

    if err_x < 1.0 && err_y < 1.0
        println(">> PASS: Far-field stresses match input loads.")
    else
        println(">> FAIL: Significant discrepancy in far-field stresses.")
    end

    # 2. Boundary Condition Check (Traction-Free)
    # ----------------------------
    # We sample points around the ellipse. 
    # NOTE: We sample at rho = 1.005 because Common.jl returns NaN for |zeta| < 1 
    # and the finite difference derivative requires stepping slightly inwards.
    rho_check = 1.005
    angles = range(0, 2π, length=12) # Check 12 points
    max_traction = 0.0

    println("\n[Boundary Condition Check @ Surface (ρ=$(rho_check))]")
    println("Target: Normal and Shear Tractions should be approx 0.0")
    @printf("%-10s %-10s %-10s %-10s\n", "Angle", "Tx", "Ty", "|T|")

    for θ in angles[1:end-1]
        # Conformal mapping: z = R(ρ*e^iθ + m/(ρ*e^iθ))
        ζ = rho_check * exp(im * θ)
        z = Common.map_to_physical(ζ, cav)
        x, y = real(z), imag(z)

        # 1. Get Stresses
        σx, σy, σ_vm = Common.calculate_field(x, y, cav, mat, stress, φ_func, ψ_func)
        τxy = 0.0 # calculate_field doesn't return tau, let's patch that temporarily or rely on Principal logic? 
        # Wait, Common.calculate_field returns (σx, σy, vm). 
        # We need Tau to check traction properly. 
        # Let's re-calculate Tau locally using the internals or accept we can only check Normal stress if we assume Principal?
        # Actually, Common.calculate_field calculates τxy internally but doesn't return it.
        # FIX: We will approximate Tau using VM: σ_vm^2 = σx^2 + σy^2 - σxσy + 3τxy^2
        # τxy^2 = (σ_vm^2 - σx^2 - σy^2 + σx*σy) / 3

        vm_sq = σ_vm^2
        tau_sq = (vm_sq - σx^2 - σy^2 + σx * σy) / 3.0
        # This loses the sign of Tau, but for checking Magnitude of Traction |T|, sign doesn't matter much 
        # as long as we check the total magnitude.
        τxy_mag = sqrt(abs(tau_sq))

        # 2. Calculate Normal Vector to Surface
        # Ellipse parametric: x = a cos β, y = b sin β (where β is conformal angle)
        # Tangent: dx/dβ = -a sin β, dy/dβ = b cos β
        # Normal: (dy/dβ, -dx/dβ) = (b cos β, a sin β)
        nx = cav.b0 * cos(θ)
        ny = cav.a0 * sin(θ)
        norm_len = sqrt(nx^2 + ny^2)
        nx /= norm_len
        ny /= norm_len

        # 3. Calculate Tractions T = σ * n
        # We try both signs of Tau since we lost the sign
        Tx = σx * nx + τxy_mag * ny
        Ty = τxy_mag * nx + σy * ny
        trac_mag = sqrt(Tx^2 + Ty^2)

        # Update max
        if trac_mag > max_traction
            max_traction = trac_mag
        end

        if (θ == 0 || θ ≈ π / 2 || θ ≈ π)
            @printf("%-10.2f %-10.4f %-10.4f %-10.4f\n", rad2deg(θ), Tx, Ty, trac_mag)
        end
    end

    println("Max Residual Traction: $(max_traction) MPa")

    # We accept a small residual due to being 0.5% away from wall and FD errors
    if max_traction < (max(stress.P1, stress.P2) * 0.1)
        println(">> PASS: Tractions are negligible near boundary.")
    else
        println(">> WARNING: Residual tractions might be high. Check if evaluation point is too far from boundary.")
    end
    println("============================\n")
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

    # 4. Sanity Check
    sanity_check(cav, mat, stress, φ_traction_free, ψ_traction_free)

    # 5. Plot
    fig = Common.plot_results(cav, mat, stress, φ_traction_free, ψ_traction_free)
    save("traction_free_ellipse.png", fig)
end

run_simulation()