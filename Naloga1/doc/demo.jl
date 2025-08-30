using LinearAlgebra, Plots, ColorSchemes
include("../src/Naloga1.jl")

"""
nariši_polinom_s_konvergenco(koeficienti; x_obseg=(-3,3), št_točk=1000, št_približkov=50)

Nariše graf polinoma, njegove ničle ter poti konvergence inverzne iteracije
za različne začetne približke.

# Argumenti
- koeficienti : koeficienti polinoma
- x_obseg     : interval za risanje grafa (privzeto (-3,3))
- št_točk     : število točk za izris polinoma (privzeto 1000)
- št_približkov : število začetnih približkov za inverzno iteracijo (privzeto 50)

# Izhod
- plt : graf
"""
function nariši_polinom_s_konvergenco(koeficienti::Vector{Float64};
        x_obseg=(-3,3), št_točk=1000, št_približkov=50)

    # companion matrika
    function companion_matrix(koef::Vector{Float64})
        n = length(koef) - 1
        C = zeros(Float64, n, n)
        for i in 1:n-1
            C[i+1, i] = 1.0
        end
        C[:, end] .= -koef[end:-1:2] ./ koef[1]
        return C
    end

    A = companion_matrix(koeficienti)
    roots = eigvals(A)

    # vrednosti polinoma
    xs = range(x_obseg[1], x_obseg[2], length=št_točk)
    ys = [sum(koeficienti[i+1] * x^(length(koeficienti)-i-1)
              for i in 0:length(koeficienti)-1) for x in xs]

    plt = plot(xs, ys, color=:black, linewidth=1, label="p(x)")

    n = length(roots)
    colors = [get(ColorSchemes.roma100, i/(n-1)) for i in 0:n-1]

    for (i, r) in enumerate(roots)
        if isreal(r)
            scatter!(plt, [real(r)], [0.0],
                     markersize=8,
                     markercolor=colors[i],
                     label="λ ≈ $(round(real(r), digits=3))")
        end
    end

    # Inverzna iteracija za začetne približke
    x0_values = range(x_obseg[1], x_obseg[2], length=št_približkov)
    line_segment_length = (x_obseg[2] - x_obseg[1]) / (3št_približkov)

    for x0 in x0_values
        λ, v = inv_lastni(A, x0)
        idx = argmin(abs.(roots .- λ))
        color = colors[idx]
        plot!(plt, [x0 - line_segment_length/2, x0 + line_segment_length/2],
              [0, 0], color=color, linewidth=3, label=false)
    end

    title!(plt, "Graf polinoma z označenimi potmi konvergence", titlefontsize=12)
    xlabel!(plt, "x")
    ylabel!(plt, "p(x)")

    return plt
end


