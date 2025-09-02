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


"""
demo_bezier(kontrolne_tocke; n=2000, nt=10_000)

Izriše Bézierjevo krivuljo in integrand ter vrne ploščino
zanke in graf obeh skupaj.

Argumenti:
- kontrolne_tocke : matrika velikosti (n+1)x2 (x,y) kontrolne točke
- n : število podintervalov za Simpsonovo pravilo (privzeto 2000)
- nt : število točk za izris krivulje (privzeto 10_000)

Vrne:
- ploščina (Float64)
- graf (Plots.Plot)
"""
function demo_bezier(kontrolne_tocke::Matrix{Float64}; n=2000, nt=10_000)
    # Parametri
    t = range(0.0, 1.0; length=nt)

    # Krivulja in odvod
    krivulja, odvodi = bezier_in_odvod(kontrolne_tocke, collect(t))
    x = krivulja[:,1]
    y = krivulja[:,2]
    dx_dt = odvodi[:,1]
    dy_dt = odvodi[:,2]

    # Integrand
    integrand = x .* dy_dt .- y .* dx_dt

    # Bézierjeva krivulja + kontrolne točke
    plt1 = plot(x, y; color=:deepskyblue1, linewidth=2, label="Bézier krivulja")
    plot!(plt1, kontrolne_tocke[:,1], kontrolne_tocke[:,2];
          seriestype=:scatter, color=:indianred1, marker=:circle, label="Kontrolne točke")
    plot!(plt1, kontrolne_tocke[:,1], kontrolne_tocke[:,2];
          color=:indianred1, linestyle=:dash, label="Povezave")
    xlabel!(plt1, "x")
    ylabel!(plt1, "y")
    title!(plt1, "Bézierjeva krivulja")
    xlims!(plt1, -1.5, 2.5)
    ylims!(plt1, -0.5, 4.5)

    # Integrand
    plt2 = plot(t, integrand; color=:indianred1, linewidth=2, label="Integrand", grid=true)
    xlabel!(plt2, "t")
    ylabel!(plt2, "Integrand")
    title!(plt2, "Integrand Bézierjeve krivulje")

    # Združeni graf
    plt = plot(plt1, plt2; layout=(1,2), size=(900,400))

    # Ploščina
    area = ploscina_bezier(kontrolne_tocke; n=n)
    println("Ploščina zanke = $(round(area, digits=12))")

    return area, plt
end

kontrolne_tocke = [
    0.0 0.0;
    1.0 1.0;
    2.0 3.0;
    1.0 4.0;
    0.0 4.0;
   -1.0 3.0;
    0.0 1.0;
    1.0 0.0
]

area, plt = demo_bezier(kontrolne_tocke; n=2000, nt=10_000)
display(plt)
