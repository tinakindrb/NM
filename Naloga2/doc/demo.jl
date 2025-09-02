using LinearAlgebra, Plots, ColorSchemes, SpecialFunctions, Statistics
include("../src/Naloga2.jl")
using .Naloga2

# Točna (referenčna) normalna CDF
F_exact(x) = 0.5 * (1 + erf(x / sqrt(2)))

# Interval vrednosti
x_vals = range(-5, 5; length=1000)

# Izračun s tvojo funkcijo in z referenco
F_aprox = [F(x) for x in x_vals]
F_ref   = [F_exact(x) for x in x_vals]

# Narišemo graf
plot(x_vals, F_ref; linestyle=:dash, linewidth=2, label="Exact F(x)")
plot!(x_vals, F_aprox; linestyle=:solid, linewidth=1, color=:red, label="Approx F(x)")
xlabel!("x")
ylabel!("F(x)")
title!("Standardna normalna CDF")

# Napaka
tolerance = 1e-10
errors = abs.(F_ref .- F_aprox)
within_tol = errors .< tolerance
percentage_ok = 100 * mean(within_tol)

println("Pravilnih izračunov: $(round(percentage_ok, digits=2)) %")
