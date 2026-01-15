using LinearAlgebra
using FFTW
using Printf
# using Revise
using Plots

# ==========================================
# 1. Structure Definitions
# ==========================================

"""
    ElasticMaterial

Holds elastic properties.
- G: Shear modulus
- ν: Poisson's ratio
- κ: Kolosov constant (3-4ν for plane strain)
"""
struct ElasticMaterial
    G::Float64
    ν::Float64
    κ::Float64
end

function ElasticMaterial(; E=15e6, ν=0.5) # Defaults from paper section 4.1
    G = E / (2 * (1 + ν))
    κ = 3 - 4 * ν # Plane strain
    return ElasticMaterial(G, ν, κ)
end

"""
    EllipticalCavity

Geometry of the cavity.
- a0, b0: Initial semi-axes
- a1, b1: Final semi-axes (deformed)
- R, m: Mapping parameters
"""
struct EllipticalCavity
    a0::Float64
    b0::Float64
    a1::Float64
    b1::Float64
    R::Float64
    m::Float64
end

function EllipticalCavity(a0, b0, a1, b1)
    R = (a0 + b0) / 2.0
    m = (a0 - b0) / (a0 + b0)
    return EllipticalCavity(a0, b0, a1, b1, R, m)
end

"""
    FarFieldStress

Stresses at infinity.
- P1, P2: Principal stresses
- Λ: Angle of P1 w.r.t x-axis
"""
struct FarFieldStress
    P1::Float64
    P2::Float64
    Λ::Float64
end

# ==========================================
# 2. Conformal Mapping & Geometry
# ==========================================

# Mapping function     z = ω(ζ, cav) [Eq. 10]
# Maps exterior of unit circle in ζ-plane to exterior of ellipse in z-plane
function ω(ζ::Complex, cav::EllipticalCavity)
    return cav.R * (ζ + cav.m / ζ)
end

# Derivative of mapping function ω'(ζ)
function dω(ζ::Complex, cav::EllipticalCavity)
    return cav.R * (1.0 - cav.m / (ζ^2))
end

# Inverse mapping z -> ζ
# Solves quadratic R*ζ^2 - z*ζ + R*m = 0
function inverse_map(z::Complex, cav::EllipticalCavity)
    # We need the root where |ζ| >= 1 (exterior)
    root1 = (z + sqrt(z^2 - 4 * cav.R^2 * cav.m)) / (2 * cav.R)
    root2 = (z - sqrt(z^2 - 4 * cav.R^2 * cav.m)) / (2 * cav.R)
    return abs(root1) >= 1.0 ? root1 : root2
end

# ==========================================
# 3. Boundary Condition Logic
# ==========================================

"""
    calculate_coefficients(cav, mat, type)

Computes Fourier coefficients A_n for the displacement boundary condition.
- type = 1: Normal displacement (Case 1)
- type = 2: Radial displacement (Case 2)
"""
function calculate_coefficients(cav::EllipticalCavity, type::Int, num_points=1024)
    ϕs = range(0, 2π, length=num_points + 1)[1:end-1]
    g_vals = ComplexF64[]

    for ϕ in ϕs
        # 1. Map Phase angle ϕ to Physical angles ϑ (Case 1) or θ (Case 2)
        # Using Eq 23 (Case 1) or Eq 30 (Case 2)
        tan_ϕ = tan(ϕ)

        # Calculate Displacement g(ϕ) = u_x + i*u_y
        u_complex = 0.0 + 0.0im

        if type == 1
            # Case 1: Normal Displacement [Source Eq 20-23]
            # Relation: tan(ϑ) = (a0^2(1-m))/(b0^2(1+m)) * tan(ϕ)
            # Careful with quadrant of atan (use atan(y, x))

            # To preserve quadrant, we expand tan terms:
            # x = R(r + m/r)cosϕ, y = R(r - m/r)sinϕ. At boundary r=1.
            x0_bound = cav.R * (1 + cav.m) * cos(ϕ)
            y0_bound = cav.R * (1 - cav.m) * sin(ϕ)

            # ϑ is angle of normal. For ellipse x = a cos t, y = b sin t, normal angle tan ϑ = a/b tan t
            # From paper Eq 23: tan ϑ = ...
            # Let's calculate ϑ directly from geometry x0, y0
            # Normal vector n = (x0/a0^2, y0/b0^2)
            nx = x0_bound / cav.a0^2
            ny = y0_bound / cav.b0^2
            ϑ = atan(ny, nx)

            # Calculate ρ(ϑ) [Eq 21 roughly - solving distance to new ellipse]
            # Simplified approach: We have (x0, y0). We need point (x1, y1) on new ellipse 
            # such that (x1-x0, y1-y0) is along normal direction ϑ.
            # This requires solving intersection. 
            # For this implementation, we approximate using the explicit ρ(ϑ) form if simpler,
            # or numerically intersect the ray from (x0, y0) along ϑ with new ellipse.

            # Ray: x = x0 + ρ cosϑ, y = y0 + ρ sinϑ
            # Subst into (x/a1)^2 + (y/b1)^2 = 1 and solve for positive ρ
            A_quad = (cos(ϑ) / cav.a1)^2 + (sin(ϑ) / cav.b1)^2
            B_quad = 2 * ((x0_bound * cos(ϑ)) / cav.a1^2 + (y0_bound * sin(ϑ)) / cav.b1^2)
            C_quad = (x0_bound / cav.a1)^2 + (y0_bound / cav.b1)^2 - 1.0

            # Quadratic formula for ρ (taking positive root closest to 0 if multiple positive, but usually one + one -)
            delta = B_quad^2 - 4 * A_quad * C_quad
            ρ = (-B_quad + sqrt(complex(delta))) / (2 * A_quad) # Take + sqrt for outward

            u_complex = ρ * exp(im * ϑ)

        elseif type == 2
            # Case 2: Radial Displacement from center [Source Eq 27-30]
            # Center angle θ
            θ = atan(y0_bound, x0_bound) # this is the physical angle

            # Ray: x = x0 + l cosθ, y = y0 + l sinθ
            # Actually, Case 2 says points move radially from center.
            # New radius r1 at angle θ on new ellipse:
            r1 = 1.0 / sqrt((cos(θ) / cav.a1)^2 + (sin(θ) / cav.b1)^2)
            r0 = sqrt(x0_bound^2 + y0_bound^2)

            l = r1 - r0
            u_complex = l * exp(im * θ)
        end

        push!(g_vals, u_complex)
    end

    # FFT to get coefficients
    # Julia FFT returns unnormalized output. Divide by N.
    # Output layout: [0, 1, 2... -2, -1] frequencies
    coeffs_raw = fft(g_vals) ./ length(g_vals)

    # We need to map these to A_n where g(ϕ) = Σ A_n σ^n
    # A_n dictionary mapping n -> value
    A_coeffs = Dict{Int,ComplexF64}()

    N = length(g_vals)
    for k in 1:N
        freq = (k <= N / 2 + 1) ? k - 1 : k - 1 - N
        A_coeffs[freq] = coeffs_raw[k]
    end

    return A_coeffs
end

# ==========================================
# 4. Complex Potentials [Source Eq 34 / 420-422]
# ==========================================

# Helper to safely get A_n (returns 0 if not present)
get_A(d::Dict{Int,ComplexF64}, n::Int) = get(d, n, 0.0im)

function φ(ζ::Complex, cav::EllipticalCavity, mat::ElasticMaterial, stress::FarFieldStress, A::Dict{Int,ComplexF64}, N_trunc::Int)
    R = cav.R
    m = cav.m
    κ = mat.κ
    P1, P2, Λ = stress.P1, stress.P2, stress.Λ

    # Term 1: -(P1+P2)/4 * R * ζ
    term1 = -((P1 + P2) / 4) * R * ζ

    # Term 2: Sum A_{-(2n-1)} / ζ^(2n-1)
    # Note: Source Eq 420 says Sum_{n=1...} A_{-(2n-1)}...
    # This implies we use negative index coefficients A_{-1}, A_{-3}, etc.
    term2 = 0.0im
    for n in 1:N_trunc
        idx = -(2 * n - 1)
        term2 += get_A(A, idx) / (ζ^(2 * n - 1))
    end
    term2 *= (2 * mat.G / κ)

    # Term 3: Far field influence in hole term
    term3_coeff = -m * (P1 + P2) / 4 + (P1 - P2) * exp(2im * Λ) / 2
    term3 = term3_coeff * (R / (κ * ζ))

    return term1 + term2 + term3
end

function ψ(ζ::Complex, cav::EllipticalCavity, mat::ElasticMaterial, stress::FarFieldStress, A::Dict{Int,ComplexF64}, N_trunc::Int)
    R = cav.R
    m = cav.m
    κ = mat.κ
    P1, P2, Λ = stress.P1, stress.P2, stress.Λ

    # Term 1: (P1-P2)e^(-2iΛ)/2 * R * ζ
    term1 = ((P1 - P2) * exp(-2im * Λ) / 2) * R * ζ

    # Term 2: -2G * Sum conjugate(A_{2n-1}) / ζ^(2n-1)
    term2 = 0.0im
    for n in 1:N_trunc
        idx = (2 * n - 1)
        term2 += conj(get_A(A, idx)) / (ζ^(2 * n - 1))
    end
    term2 *= -2 * mat.G

    # Term 3
    term3 = R * ((P1 + P2) / 4) * (κ / ζ - ζ * (1 + m^2) / (ζ^2 - m))

    # Term 4: Complex interaction term
    # -(ζ * (1 + mζ^2)/(ζ^2 - m)) * { ... }
    brackets = 0.0im

    # Part A inside brackets: Sum (2n-1) * A_{-(2n-1)} / ζ^(2n)
    sum_part = 0.0im
    for n in 1:N_trunc
        idx = -(2 * n - 1)
        sum_part += (2 * n - 1) * get_A(A, idx) / (ζ^(2 * n))
    end
    sum_part *= (2 * mat.G / κ)

    # Part B inside brackets
    partB = (m * (P1 + P2) / 4 - (P1 - P2) * exp(2im * Λ) / 2) * (R / (κ * ζ^2))

    term4 = -(ζ * (1 + m * ζ^2) / (ζ^2 - m)) * (sum_part + partB)

    return term1 + term2 + term3 + term4
end

# Derivatives for Stress Calculation Φ(ζ) = φ'(ζ)/ω'(ζ)
# We calculate φ'(ζ) analytically from the series above
function dφ(ζ::Complex, cav::EllipticalCavity, mat::ElasticMaterial, stress::FarFieldStress, A::Dict{Int,ComplexF64}, N_trunc::Int)
    R = cav.R
    m = cav.m
    κ = mat.κ
    P1, P2, Λ = stress.P1, stress.P2, stress.Λ

    # 1. Derivative of linear term
    d_term1 = -((P1 + P2) / 4) * R

    # 2. Derivative of sum term A_{-(2n-1)} * ζ^(-(2n-1))
    # d/dζ (ζ^-k) = -k * ζ^(-k-1)
    d_term2 = 0.0im
    for n in 1:N_trunc
        idx = -(2 * n - 1)
        pow = 2 * n - 1
        d_term2 += get_A(A, idx) * (-pow) / (ζ^(pow + 1))
    end
    d_term2 *= (2 * mat.G / κ)

    # 3. Derivative of 1/ζ term
    term3_coeff = -m * (P1 + P2) / 4 + (P1 - P2) * exp(2im * Λ) / 2
    d_term3 = term3_coeff * (R / κ) * (-1 / ζ^2)

    return d_term1 + d_term2 + d_term3
end

function dψ(ζ::Complex, cav::EllipticalCavity, mat::ElasticMaterial, stress::FarFieldStress, A::Dict{Int,ComplexF64}, N_trunc::Int)
    # For stress, we need Ψ(ζ) = ψ'(ζ)/ω'(ζ)
    # We can approximate derivative numerically or derive analytically. 
    # Given the complexity of term 4 in ψ, a finite difference is often robust enough for engineering scripts
    # unless exact analytical precision is required.
    h = 1e-6
    return (ψ(ζ + h, cav, mat, stress, A, N_trunc) - ψ(ζ - h, cav, mat, stress, A, N_trunc)) / (2 * h)
end

# ==========================================
# 5. Solver
# ==========================================

function solve_at_point(x::Float64, y::Float64, cav::EllipticalCavity, mat::ElasticMaterial, stress::FarFieldStress, A_coeffs::Dict, N_trunc=50)
    z = x + im * y
    ζ = inverse_map(z, cav)

    # Check if inside cavity (approximate)
    if abs(ζ) < 1.0
        return NaN, NaN, NaN, NaN, NaN
    end

    # Compute Potentials and Derivatives
    phi_val = φ(ζ, cav, mat, stress, A_coeffs, N_trunc)
    psi_val = ψ(ζ, cav, mat, stress, A_coeffs, N_trunc)

    dphi_val = dφ(ζ, cav, mat, stress, A_coeffs, N_trunc)
    dpsi_val = dψ(ζ, cav, mat, stress, A_coeffs, N_trunc) # Numerical derivative

    omega_prime = dω(ζ, cav)

    Phi = dphi_val / omega_prime
    Psi = dpsi_val / omega_prime

    # Stresses [Source Eq 179-180]
    # σx + σy = 4 Re[Φ(ζ)]
    sum_sigma = 4 * real(Phi)

    # σy - σx + 2iτxy = 2 [ conj(ω(ζ))/ω'(ζ) * Φ'(ζ) + Ψ(ζ) ]
    # Wait, the formula in paper Eq 8 (Source 128) uses z_bar.
    # In mapped plane (Source Eq 180): 2 [ (conj(ω(ζ))/ω'(ζ)) * Φ'(ζ) + Ψ(ζ) ]
    # Note: Φ'(ζ) here means d/dz Φ? No, usually formulas are in ζ.
    # Standard Muskhelishvili in ζ:
    # σy - σx + 2iτxy = 2 * [ (conj(ω(ζ))/ω'(ζ)) * (dΦ/dζ / ω'(ζ))??? ] -> No.
    # Let's check Source 180 carefully:
    # " ... = 2 [ (conj(ω(ζ))/ω'(ζ)) * Φ'(ζ) + Ψ(ζ) ] "
    # where Φ(ζ) = φ'(ζ)/ω'(ζ) and Ψ(ζ) = ψ'(ζ)/ω'(ζ).
    # The term Φ'(ζ) in the brackets usually refers to d(Φ)/dζ.

    # Calculate dΦ/dζ numerically
    h = 1e-6
    dphi_plus = dφ(ζ + h, cav, mat, stress, A_coeffs, N_trunc)
    omega_prime_plus = dω(ζ + h, cav)
    Phi_plus = dphi_plus / omega_prime_plus

    dphi_minus = dφ(ζ - h, cav, mat, stress, A_coeffs, N_trunc)
    omega_prime_minus = dω(ζ - h, cav)
    Phi_minus = dphi_minus / omega_prime_minus

    dPhi_dzeta = (Phi_plus - Phi_minus) / (2 * h)

    # Now assemble
    term_diff = 2 * ((conj(ω(ζ, cav)) / omega_prime) * dPhi_dzeta + Psi)

    sigma_y_minus_x = real(term_diff)
    tau_xy_2 = imag(term_diff)

    sigma_x = 0.5 * (sum_sigma - sigma_y_minus_x)
    sigma_y = 0.5 * (sum_sigma + sigma_y_minus_x)
    tau_xy = 0.5 * tau_xy_2

    # Displacements [Source Eq 180 / 9]
    # 2G(ux + iuy) = κ φ(ζ) - (ω(ζ)/conj(ω'(ζ))) * conj(φ'(ζ)) - conj(ψ(ζ))
    # Note: Eq 9 uses z_bar. Eq 180 uses mapped version.
    # Source 180: ... = κ φ(ζ) - (ω(ζ)/conj(ω'(ζ))) * conj(φ'(ζ)) - conj(ψ(ζ))

    disp_rhs = mat.κ * phi_val - (ω(ζ, cav) / conj(omega_prime)) * conj(dphi_val) - conj(psi_val)
    disp = disp_rhs / (2 * mat.G)

    return sigma_x, sigma_y, tau_xy, real(disp), imag(disp)
end

# ==========================================
# 6. Example Run (Figure 8 Params)
# ==========================================

function run_example()
    # Parameters from Figure 8 in the paper [Source 559-561]
    # a0=5, b0=3, a1=5.1, b1=3.05 (units mm, but consistent units work)
    # P1=500 kPa (x-axis), P2=250 kPa (y-axis)
    # G=20 MPa = 20000 kPa, nu=0.4

    a0, b0 = 5.0, 3.0
    a1, b1 = 5.1, 3.05
    G_mod = 20000.0 # kPa
    nu = 0.4

    P1 = 500.0
    P2 = 250.0
    Lambda = 0.0 # P1 along x-axis

    cavity = EllipticalCavity(a0, b0, a1, b1)
    material = ElasticMaterial(G_mod, nu, 3 - 4 * nu)
    stress = FarFieldStress(P1, P2, Lambda)

    # Calculate Coefficients for Case 1 (Normal Displacement)
    println("Calculating Fourier Coefficients for Case 1...")
    coeffs = calculate_coefficients(cavity, 1) # 1 for Normal

    # --- 1D Plotting along X-axis (as before) ---
    println("Calculating Stresses along X-axis for 1D plot...")
    xs_1d = range(a0, 3 * a0, length=100)
    sig_y_vals = Float64[]
    u_x_vals = Float64[]
    
    for x in xs_1d
        sx, sy, txy, ux, uy = solve_at_point(x, 0.0, cavity, material, stress, coeffs)
        push!(sig_y_vals, sy)
        push!(u_x_vals, ux)
    end

    p_stress_1d = plot(xs_1d ./ a0, sig_y_vals ./ P1,
        label="σy / P1",
        xlabel="x / a0",
        ylabel="Stress Concentration",
        title="Stress along x-axis",
        lw=2)
    hline!([stress.P2 / stress.P1], linestyle=:dash, label="P2/P1 (Far-field)")

    p_disp_1d = plot(xs_1d ./ a0, u_x_vals ./ a0,
        label="ux / a0",
        xlabel="x / a0",
        ylabel="Normalized Displacement",
        title="Displacement along x-axis",
        lw=2,
        color=:red)

    plot_1d = plot(p_stress_1d, p_disp_1d, layout=(2, 1), legend=:best)
    display(plot_1d)
    println("1D plots displayed. Now generating 2D field plots...")

    # --- 2D Field Plotting ---
    grid_res = 100
    lim = 2.5 * a0
    xs_2d = range(-lim, lim, length=grid_res)
    ys_2d = range(-lim, lim, length=grid_res)

    sigma_vm = zeros(grid_res, grid_res)
    disp_mag = zeros(grid_res, grid_res)

    for (i, y) in enumerate(ys_2d), (j, x) in enumerate(xs_2d)
        sx, sy, txy, ux, uy = solve_at_point(x, y, cavity, material, stress, coeffs)
        
        if isnan(sx)
            sigma_vm[i, j] = NaN
            disp_mag[i, j] = NaN
        else
            # Von Mises Stress
            sigma_vm[i, j] = sqrt(sx^2 - sx*sy + sy^2 + 3*txy^2)
            # Displacement Magnitude
            disp_mag[i, j] = sqrt(ux^2 + uy^2)
        end
    end

    # Function to draw the ellipse boundary
    function ellipse_shape(cav, n=100)
        t = range(0, 2π, length=n)
        x = cav.a0 .* cos.(t)
        y = cav.b0 .* sin.(t)
        return Shape(x, y)
    end

    p_stress_2d = heatmap(xs_2d, ys_2d, sigma_vm,
        aspect_ratio=:equal,
        c=:viridis,
        title="Von Mises Stress Field",
        xlabel="x",
        ylabel="y")
    plot!(p_stress_2d, ellipse_shape(cavity), fillalpha=0, lw=2, linecolor=:white, label="Cavity")

    p_disp_2d = heatmap(xs_2d, ys_2d, disp_mag,
        aspect_ratio=:equal,
        c=:inferno,
        title="Displacement Magnitude Field",
        xlabel="x",
        ylabel="y")
    plot!(p_disp_2d, ellipse_shape(cavity), fillalpha=0, lw=2, linecolor=:white, label="Cavity")

    plot_2d = plot(p_stress_2d, p_disp_2d, layout=(1, 2), size=(1200, 500))
    display(plot_2d)
    println("2D field plots displayed.")
run_example()