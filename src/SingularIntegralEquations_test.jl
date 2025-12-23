# Using SingularIntegralEquations for singular integrals
using ApproxFun, SingularIntegralEquations

# Define the domain and function space
d = ChebyshevInterval()  # [-1, 1] interval
S = JacobiWeight(0.5, 0.5, Chebyshev(d))  # Weighted Chebyshev space

# Define a test function, e.g., f(t) = exp(t)
x = Fun(d)
f = Fun(exp, S)

# Define a point 'a' inside the interval at which to evaluate the integral
a = 0.1

# Compute the Cauchy Principal Value (Hilbert Transform):
# H(f)(a) = (1/π) * ⨎_{-1}^{1} f(t) / (t - a) dt
H = hilbert(f, a)
println("Hilbert transform of exp(x) at x=$a: ", H)

# Compute an integral with a logarithmic singularity:
# L(f)(a) = ∫_{-1}^{1} f(t) * log|t - a| dt
L = logkernel(f, a)
println("Log kernel integral of exp(x) at x=$a: ", L)

# Compute the Stieltjes integral:
# S(f)(a) = (1/π) * ∫_{-1}^{1} f(t) / (a - t) dt
St = stieltjes(f, Complex(a))
println("Stieltjes integral of exp(x) at x=$a: ", St)