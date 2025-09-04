module Naloga2

export simpsonovo_pravilo, F, binomski, bezier, bezier_in_odvod, ploscina_bezier

"""
simpsonovo_pravilo(f, a, b, n)

Izračuna približek določenega integrala funkcije f na intervalu [a,b] z 
uporabo sestavljenega Simpsonovega pravila.

Vhod:
- f : funkcija, ki jo integriramo
- a : začetek intervala integracije
- b : konec intervala integracije
- n : število podintervalov (mora biti sodo)

Izhod: 
- približek integrala
"""
function simpsonovo_pravilo(f, a::Float64, b::Float64, n::Int)
    if n % 2 != 0
        error("Število intervalov n mora biti sodo.")
    end
    h = (b - a) / n
    S = f(a) + f(b)
    for i in 1:2:(n-1)
        S += 4 * f(a + i*h)
    end
    for i in 2:2:(n-2)
        S += 2 * f(a + i*h)
    end
    return S * h / 3
end


"""
F(x; n=2000)

Porazdelitvena funkcija standardne normalne porazdelitve N(0,1)

Metoda:
- Za x > 0: integral na intervalu [0, x] z uporabo Simpsonovega pravila.
- Za x < 0: uporabi simetrijo F(x) = 1 - F(-x).
- Za zelo velike |x|: uporabi asimptotsko aproksimacijo repa porazdelitve.
- Za x = 0: F(0) = 0.5.

Vhod:
- x : vrednost, pri kateri računamo porazdelitveno funkcijo
- n: število podintervalov za Simpsonovo pravilo (privzeto 2000)

Izhod:
- približek F(x)
"""
function F(x::Float64; n=2000)
    pdf(t) = exp(-t^2/2) / sqrt(2π)  # gostota normalne porazdelitve
    
    if x < -8
        return 0.0
    elseif x > 8
        return 1.0 - pdf(x) / (x + 1/x)
    elseif x == 0.0
        return 0.5
    elseif x > 0
        return 0.5 + simpsonovo_pravilo(pdf, 0.0, x, n)
    else
        return 1.0 - F(-x; n=n)
    end
end


"""
binomski(n, k)

Izračuna binomski koeficient "n nad k".
"""
function binomski(n::Int, k::Int)
    if k < 0 || k > n
        return 0
    end
    k = min(k, n-k)
    c = 1
    for i in 1:k
        c = c * (n - (k - i)) ÷ i
    end
    return c
end


"""
bezier(kontrolne_tocke, t_vrednosti)

Izračuna točke Bézierjeve krivulje za dane kontrolne točke in parametre.

Vhod:
- kontrolne_tocke : matrika velikosti (n+1)x2, vsaka vrstica je točka (x,y).
- t_vrednosti : vektor vrednosti parametra t ∈ [0,1].

Izhod:
- matrika velikosti length(t_vrednosti)x2 z izračunanimi točkami na krivulji.
"""
function bezier(kontrolne_tocke::Matrix{Float64}, t_vrednosti::Vector{Float64})
    n = size(kontrolne_tocke, 1) - 1
    krivulja = zeros(length(t_vrednosti), 2)

    for i in 0:n
        binomski_koef = binomski(n, i)
        bernsteinov_polinom = binomski_koef .* (t_vrednosti.^i) .* ((1 .- t_vrednosti).^(n-i))
        krivulja .+= bernsteinov_polinom .* kontrolne_tocke[i+1, :]'
    end
    return krivulja
end


"""
bezier_in_odvod(kontrolne_tocke, t_vrednosti)

Izračuna točke Bézierjeve krivulje in njen odvod
za dane kontrolne točke in parametre.

Vhod:
- kontrolne_tocke : matrika velikosti (n+1)×2, vsaka vrstica je točka (x,y).
- t_vrednosti : vektor vrednosti parametra t ∈ [0,1].

Izhod:
- matriko točk na krivulji
- matriko točk odvoda (x'(t), y'(t))
"""
function bezier_in_odvod(kontrolne_tocke::Matrix{Float64}, t_vrednosti::Vector{Float64})
    n = size(kontrolne_tocke, 1) - 1

    # krivulja stopnje n
    krivulja = bezier(kontrolne_tocke, t_vrednosti)

    # kontrolne točke za odvod (razlike med zaporednimi)
    odvodne_tocke = [n .* (kontrolne_tocke[i+1,:] .- kontrolne_tocke[i,:]) for i in 1:n]
    odvodne_tocke = permutedims(hcat(odvodne_tocke...))  # velikost n × 2 (prava matrika)

    # odvod je Bézierjeva krivulja stopnje n-1
    odvodi = bezier(odvodne_tocke, t_vrednosti)

    return krivulja, odvodi
end


"""
ploscina_bezier(kontrolne_tocke; n=1000)

Izračuna ploščino, ki jo omejuje Bézierjeva krivulja,
z uporabo formule 1/2 ∫ (x(t)*y'(t) - x'(t)*y(t)) dt.

Vhod:
- kontrolne_tocke : matrika velikosti (m)x2, kontrolne točke (x,y).
- n : število podintervalov za Simpsonovo pravilo (privzeto 1000).

Izhod:
- ploščina.
"""
function ploscina_bezier(kontrolne_tocke::Matrix{Float64}; n::Int=1000)
    t_vrednosti = range(0.0, 1.0; length=n+1)

    krivulja, odvodi = bezier_in_odvod(kontrolne_tocke, collect(t_vrednosti))
    x = krivulja[:,1]
    y = krivulja[:,2]
    dx_dt = odvodi[:,1]
    dy_dt = odvodi[:,2]

    integrand(t) = begin
        i = Int(round(t*n)) + 1
        return x[i]*dy_dt[i] - y[i]*dx_dt[i]
    end

    integral = simpsonovo_pravilo(integrand, 0.0, 1.0, n)
    return 0.5 * integral
end

end # module Naloga2
