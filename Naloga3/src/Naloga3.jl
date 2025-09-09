module Naloga3

using LinearAlgebra

export dopri5, sistem_nihalo, energija_nihalo, resi_nihalo

"""
dopri5(f, tspan, y0, epsilon, sigma)

Reši začetni problem za navadne diferencialne enačbe z metodo Dormand-Prince 
5. reda (DOPRI5) z adaptivnim korakom in oceno lokalne napake.

Vhod:
- f : funkcija desne strani sistema (dy/dt = f(t, y)), kjer je y vektor
- tspan : interval reševanja [t0, tf]
- y0 : začetni pogoj (vektor)
- epsilon : toleranca za lokalno napako
- sigma : varnostni faktor (0 < sigma < 1) pri prilagoditvi koraka

Izhod:
- T : vektor časovnih točk, kjer je rešitev izračunana
- Y : matrika rešitev, kjer vrstica ustreza času v T
"""
function dopri5(f, tspan, y0, epsilon, sigma)
    t0, tf = tspan
    t = t0
    y = copy(y0)
    h = (tf - t0) / 100 # začetni korak
    T = [t]
    Y = [copy(y)]

    while t < tf
        k1 = h * f(t, y)
        k2 = h * f(t + h/5, y + k1/5)
        k3 = h * f(t + 3h/10, y + 3k1/40 + 9k2/40)
        k4 = h * f(t + 4h/5, y + 44k1/45 - 56k2/15 + 32k3/9)
        k5 = h * f(t + 8h/9, y + 19372k1/6561 - 25360k2/2187 + 64448k3/6561 - 212k4/729)
        k6 = h * f(t + h, y + 9017k1/3168 - 355k2/33 + 46732k3/5247 + 49k4/176 - 5103k5/18656)
        k7 = h * f(t + h, y + 35k1/384 + 500k3/1113 + 125k4/192 - 2187k5/6784 + 11k6/84)

        # Lokalna napaka (4. red)
        le = 71k1/57600 - 71k3/16695 + 71k4/1920 - 17253k5/339200 + 22k6/525 - k7/40
        nle = norm(le)

        if nle < epsilon * h
            # Sprejmemo korak (5. red)
            y += 35k1/384 + 500k3/1113 + 125k4/192 - 2187k5/6784 + 11k6/84
            t += h
            push!(T, t)
            push!(Y, copy(y))

            # Prilagoditev koraka
            h = h * sigma * (epsilon * h / max(nle, 1e-16))^(1/5)
            if t + h > tf
                h = tf - t
            end
        else
            h /= 2
        end
    end

    return T, hcat(Y...)'
end

"""
sistem_nihalo(t, y)

Definira sistem diferencialnih enačb za matematično nihalo:
- kot'(t) = kotna_hitrost(t)
- kotna_hitrost'(t) = -(g/l) * sin(kot(t))

Vhod:
- t : čas
- y : vektor [kot, kotna_hitrost]

Izhod:
- odvod vektorja [kot, kotna_hitrost]
"""
function sistem_nihalo(t, y)
    kot, kotna_hitrost = y
    g, l = 9.81, 1.0
    return [kotna_hitrost, -(g/l) * sin(kot)]
end

"""
energija_nihalo(y)

Izračuna mehansko energijo matematičnega nihala za dano stanje.

Energija je sestavljena iz kinetične in potencialne energije:
E = 0.5 * kotna_hitrost^2 + (g/l) * (1 - cos(kot))

Vhod:
- y : vektor [kot, kotna_hitrost]

Izhod:
- energija : skalarna vrednost energije
"""
function energija_nihalo(y)
    kot, kotna_hitrost = y
    g, l = 9.81, 1.0
    return 0.5 * kotna_hitrost^2 + (g/l) * (1 - cos(kot))
end

"""
resi_nihalo(zacetni_kot, zacetna_kotna_hitrost, cas_konec; epsilon=1e-6, sigma=0.9)

Reši matematično nihalo z metodo DOPRI5 za dane začetne pogoje in časovni interval.

Vhod:
- zacetni_kot : začetni odmik (v radianih)
- zacetna_kotna_hitrost : začetna kotna hitrost
- cas_konec : čas integracije
- epsilon : toleranca napake (privzeto 1e-6)
- sigma : varnostni faktor (privzeto 0.9)

Izhod:
- cas : vektor časov
- resitev : matrika [kot, kotna_hitrost] za vsak čas
"""
function resi_nihalo(zacetni_kot, zacetna_kotna_hitrost, cas_konec; epsilon=1e-6, sigma=0.9)
    zacetni_pogoj = [zacetni_kot, zacetna_kotna_hitrost]
    casovni_interval = (0.0, cas_konec)
    return dopri5(sistem_nihalo, casovni_interval, zacetni_pogoj, epsilon, sigma)
end

end # module Naloga3
