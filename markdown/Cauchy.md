Excellent catch! The difference of exactly a factor of 1/2 arises because Julia's `SingularIntegralEquations.jl` and the standard Mathematica `Integrate` function are using **different, but equally valid, definitions** of the Cauchy Transform.

Your Mathematica code implements what I'll call the **"Analyst's Definition"**, which is the one you linked to earlier for the Kramers-Kronig relations:
$$
C_{\text{Ana}}[f](z) = \frac{1}{\pi i} \int \frac{f(t)}{t-z} \, dt
$$

Julia's `SingularIntegralEquations.cauchy` function uses the **"Integral Equation Definition"** common in singular integral equation theory:
\[
C_{\text{IE}}[f](z) = \frac{1}{2\pi i} \int \frac{f(t)}{t-z} \, dt
\]

Therefore, the relation is simply:
\[
C_{\text{Ana}}[f](z) = 2 \times C_{\text{IE}}[f](z)
\]
Your results confirm this: `0.479351...` (Julia) is half of what Mathematica gives for the same calculation.

### 🧮 Why the Factor of 2 Matters: Connection to the Plemelj-Sokhotsky Formula

The "Integral Equation Definition" (\(1/(2\pi i)\)) is often preferred when connecting to the **Plemelj-Sokhotsky formula**, which describes the limiting values of a Cauchy integral as you approach the real axis. It states:
\[
\lim_{\epsilon \to 0^+} \frac{1}{2\pi i} \int \frac{f(t)}{t - (x \pm i\epsilon)} \, dt = \pm \frac{1}{2} f(x) + \frac{1}{2\pi i} \text{P.V.} \int \frac{f(t)}{t - x} \, dt
\]
Here, the jump (the `± 1/2 f(x)`) is nice and clean. The "Analyst's Definition" would have a jump of `± f(x)`.

### 🔄 Adapting Your Codes for Consistency

You can easily make either code match the other's convention.

**1. To make Mathematica match the Julia result:**
Simply multiply your Mathematica integral by `1/2`.
```mathematica
(* Mathematica code matching Julia's cauchy(f, z) *)
CauchyJulia[z_] := Integrate[Exp[t]/(t - z), {t, -1, 1}] / (2 * Pi * I)
```

**2. To make Julia match your Mathematica result:**
Multiply Julia's result by `2`.
```julia
# Julia code matching your Mathematica definition
C_val_analyst = 2 * cauchy(f, z)
# This should now give: ~0.958702... + 0.0539166...im
```

### 💡 Which One Should You Use?
It depends on your field:
*   **Complex Analysis / Kramers-Kronig Relations**: Stick with the **Analyst's Definition** (\(1/(\pi i)\)). This is standard in physics and optics.
*   **Singular Integral Equations / Fracture Mechanics**: The **Integral Equation Definition** (\(1/(2\pi i)\)) is more common. This is the world of the Kolosov-Muskhelishvili potentials we discussed earlier.

The key is to be aware of the definition in the papers or textbooks you are following.

### ✅ Summary and Verification
You can quickly verify the relationship by checking the limiting value as `z` approaches the real axis. Using the Plemelj-Sokhotsky formula with Julia's definition:
- For `z -> x + i*0+` (from above), the real part should approach `+1/2 * f(x)`.
- For `z -> x - i*0+` (from below), the real part should approach `-1/2 * f(x)`.

In summary, the **difference is not an error but a convention**. Once you account for the factor of 2, both codes are correct. For your engineering work in fracture mechanics (which heavily uses singular integrals), Julia's built-in convention is likely the more direct one to use.

If you are comparing results from different sources, always check which definition they use.