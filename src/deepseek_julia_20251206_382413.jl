# Import the necessary packages
using Pkg
Pkg.add(["ApproxFun", "SingularIntegralEquations"])
using ApproxFun, SingularIntegralEquations

# Define the interval and the function space (Chebyshev polynomials on [-1, 1])
x = Fun()
# Define a test function, e.g., f(t) = exp(t)
f = exp(x)

# Define a point 'a' inside the interval at which to evaluate the integral
a = 0.1

# Compute the Cauchy Principal Value (Hilbert Transform):
# H(f)(a) = (1/π) * ⨎_{-1}^{1} f(t) / (t - a) dt
H = hilbert(f, a)
println("Hilbert transform of exp(x) at x=$a: ", H)

# Compute an integral with a logarithmic singularity (appears in fracture mechanics):
# L(f)(a) = (1/π) * ∫_{-1}^{1} f(t) * log|t - a| dt
L = logkernel(f, a)
println("Log kernel integral of exp(x) at x=$a: ", L)