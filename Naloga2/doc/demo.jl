using LinearAlgebra, Plots, ColorSchemes, SpecialFunctions, Statistics
include("../src/Naloga2.jl")
using .Naloga2

"""
    demo_F(; tolerance=1e-10, npoints=1000)

Primerjava tvoje implementacije funkcije F z referenčno vrednostjo
(izračunano z `erf` iz SpecialFunctions). Nariše graf in vrne graf + 
delež točk, ki so znotraj tolerance.

Argumenti:
- tolerance::Float64 : meja za napako (privzeto 1e-10)
- npoints::Int : število točk v intervalu [-5,5] (privzeto 1000)

Vrne:
- (plt, percentage_ok) kjer je `plt` graf in `percentage_ok` procent pravilnih točk.
"""
function demo_F(; tolerance=1e-10, npoints=1000)
    # Točna (referenčna) normalna CDF
    F_exact(x) = 0.5 * (1 + erf(x / sqrt(2)))

    # Interval vrednosti
    x_vals = range(-5, 5; length=npoints)

    # Izračun s tvojo funkcijo in z referenco
    F_aprox = [F(x) for x in x_vals]
    F_ref   = [F_exact(x) for x in x_vals]

    # Graf
    plt = plot(x_vals, F_ref; linestyle=:dash, linewidth=2, label="Exact F(x)")
    plot!(plt, x_vals, F_aprox; linestyle=:solid, linewidth=1, color=:red, label="Approx F(x)")
    xlabel!("x")
    ylabel!("F(x)")
    title!("Standardna normalna CDF")

    # Napaka
    errors = abs.(F_ref .- F_aprox)
    within_tol = errors .< tolerance
    percentage_ok = 100 * mean(within_tol)

    println("Pravilnih izračunov: $(round(percentage_ok, digits=2)) %")

    return plt, percentage_ok
end

plt, pct = demo_F()
display(plt)
println("Ujemanje: $pct %")
