module Naloga2

export simpsonovo_pravilo, F

"""
simpsonovo_pravilo(f, a, b, n)

Izračuna približek določenega integrala funkcijef na intervalu [a,b] z 
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
- x : vrednost, pri kateri računamo CDF
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


end # module Naloga2
