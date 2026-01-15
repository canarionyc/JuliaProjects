# Common.jl
module Common

using CairoMakie, Makie
using LinearAlgebra
using Printf

# ====================
# Data Structures
# ====================

struct EllipticalCavity
    a0::Float64
    b0::Float64
    R::Float64
    m::Float64
end

function EllipticalCavity(a0, b0)
    R = (a0 + b0) / 2.0
    m = (a0 - b0) / (a0 + b0)
    return EllipticalCavity(a0, b0, R, m)
end

struct ElasticMaterial
    G::Float64
    ν::Float64
    κ::Float64
end

function ElasticMaterial(; E=15e6, ν=0.5)
    G = E / (2 * (1 + ν))
    κ = 3 - 4 * ν # Plane strain
    return ElasticMaterial(G, ν, κ)
end

struct FarFieldStress
    P1::Float64
    P2::Float64
    Λ::Float64
end

# ====================
# Geometry / Mapping
# ====================

# z = ω(ζ) = R(ζ + m/ζ)
function map_to_physical(ζ::Complex, cav::EllipticalCavity)
    return cav.R * (ζ + cav.m / ζ)
end

# Derivative ω'(ζ) = R(1 - m/ζ²)
function d_map(ζ::Complex, cav::EllipticalCavity)
    return cav.R * (1.0 - cav.m / (ζ^2))
end

# Inverse z -> ζ
function map_to_phase(z::Complex, cav::EllipticalCavity)
    # Roots of R*ζ² - z*ζ + R*m = 0
    # We take the root with |ζ| ≥ 1 (exterior)
    disc = sqrt(z^2 - 4 * cav.R^2 * cav.m)
    r1 = (z + disc) / (2 * cav.R)
    r2 = (z - disc) / (2 * cav.R)
    return abs(r1) >= 1.0 ? r1 : r2
end

# ====================
# Solver Helper
# ====================

# Generic stress calculator using potentials provided by the specific problem file
# φ_func and ψ_func must signature: (ζ, cav, mat, stress) -> Complex
function calculate_field(x, y, cav, mat, stress, φ_func, ψ_func)
    z = x + im*y
    ζ = map_to_phase(z, cav)

    # Mask inside cavity (approximately |ζ| < 1)
    if abs(ζ) < 1.0
        return NaN, NaN, NaN
    end

    # Evaluate Potentials
    # Note: derivatives dφ/dζ are usually needed. 
    # For robustness in these exercises, we can use Finite Difference for derivatives
    # unless analytical dφ is provided.
    h = 1e-6
    val_φ = φ_func(ζ, cav, mat, stress)
    der_φ = (φ_func(ζ+h, cav, mat, stress) - φ_func(ζ-h, cav, mat, stress)) / (2h)
    der_ψ = (ψ_func(ζ+h, cav, mat, stress) - ψ_func(ζ-h, cav, mat, stress)) / (2h)
    
    ω_prime = d_map(ζ, cav)
    
    # Kolmogorov-Muskhelishvili Stresses (in ζ plane)
    # Φ(ζ) = φ'(ζ) / ω'(ζ)
    # Ψ(ζ) = ψ'(ζ) / ω'(ζ)
    Phi = der_φ / ω_prime
    Psi = der_ψ / ω_prime

    # σx + σy = 4 Re[Φ]
    # σy - σx + 2iτxy = 2 [ conj(ω(ζ))/ω'(ζ) * Φ'(ζ) + Ψ(ζ) ]
    
    # We need Φ'(ζ) (derivative wrt z? No, usually wrt z in the formula).
    # Standard formula: σy - σx + 2iτxy = 2 [ z_bar * Φ'(z) + Ψ(z) ]
    # Transformed: 2 [ (conj(ω)/ω') * d(Φ)/dζ + Ψ ] ??? 
    # Let's use the robust form:
    # Term2 = 2 * ( conj(map_to_physical(ζ, cav)) * (der_Phi_z) + Psi )
    # where der_Phi_z = dΦ/dζ * (1/ω')
    
    # Calculate dΦ/dζ numerically
    Phi_plus  = ((φ_func(ζ+h, cav, mat, stress) - φ_func(ζ, cav, mat, stress))/h) / d_map(ζ+h, cav)
    Phi_minus = ((φ_func(ζ, cav, mat, stress) - φ_func(ζ-h, cav, mat, stress))/h) / d_map(ζ-h, cav)
    dPhi_dzeta = (Phi_plus - Phi_minus) / h # approximate central diff of Phi
    
    der_Phi_z = dPhi_dzeta / ω_prime

    sum_s = 4 * real(Phi)
    diff_s = 2 * ( (conj(map_to_physical(ζ, cav)) / ω_prime) * dPhi_dzeta + Psi ) # Check formula consistency

    # Note: The term often is conj(ω)/ω' * dΦ/dζ. 
    # Let's stick to the form used in the previous solution for consistency.
    
    σx = real(sum_s - diff_s) / 2.0
    σy = real(sum_s + diff_s) / 2.0
    τxy = imag(diff_s) / 2.0
    
    # Von Mises
    σ_vm = sqrt(σx^2 + σy^2 - σx*σy + 3*τxy^2)
    
    return σx, σy, σ_vm
end

# ====================
# Plotting
# ====================

function plot_results(cav, mat, stress, φ_func, ψ_func; limit=3.0, res=200)
    xs = range(-limit*cav.a0, limit*cav.a0, length=res)
    ys = range(-limit*cav.a0, limit*cav.a0, length=res)

    σx_grid = zeros(res, res)
    σy_grid = zeros(res, res)
    vm_grid = zeros(res, res)

    for (j, y) in enumerate(ys)
        for (i, x) in enumerate(xs)
            sx, sy, vm = calculate_field(x, y, cav, mat, stress, φ_func, ψ_func)
            σx_grid[j, i] = sx
            σy_grid[j, i] = sy
            vm_grid[j, i] = vm
        end
    end

    # Determine shared color limits for Sigma X and Sigma Y
    valid_stresses = filter(!isnan, [σx_grid..., σy_grid...])
    min_stress = minimum(valid_stresses)
    max_stress = maximum(valid_stresses)
    clims_sigma = (min_stress, max_stress)
    
    clims_vm = (0, maximum(filter(!isnan, vm_grid))*0.8)

    fig = Figure(size = (800, 900))

    ax1 = Axis(fig[1, 1], aspect = DataAspect(), title = "Sigma X")
    hm1 = heatmap!(ax1, xs, ys, σx_grid, colormap = :viridis, colorrange = clims_sigma)
    
    ax2 = Axis(fig[1, 2], aspect = DataAspect(), title = "Sigma Y")
    heatmap!(ax2, xs, ys, σy_grid, colormap = :viridis, colorrange = clims_sigma)

    ax3 = Axis(fig[2, 1], aspect = DataAspect(), title = "Von Mises")
    hm3 = heatmap!(ax3, xs, ys, vm_grid, colormap = :inferno, colorrange = clims_vm)

    # Hide the empty axis
    ax4 = Axis(fig[2, 2], aspect = DataAspect())
    hidedecorations!(ax4)
    hidespines!(ax4)

    Colorbar(fig[3, 1:2], hm1, label = "Normal Stress (σx, σy)", vertical = false, flipaxis = false)
    Colorbar(fig[4, 1:2], hm3, label = "Von Mises Stress", vertical = false, flipaxis = false)


    # Draw the ellipse outline
    θ = range(0, 2π, length=100)
    el_x = cav.a0 .* cos.(θ)
    el_y = cav.b0 .* sin.(θ)
    
    lines!(ax1, el_x, el_y, color=:white, linewidth=2)
    lines!(ax2, el_x, el_y, color=:white, linewidth=2)
    lines!(ax3, el_x, el_y, color=:white, linewidth=2)

    return fig
end


end # module
