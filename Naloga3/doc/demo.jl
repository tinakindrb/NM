using Plots
include("../src/Naloga3.jl")
using .Naloga3

# ------------------------------
# 1. Linearno (harmonično) nihalo
# ------------------------------
function linearno_nihalo(cas, kot0, kotna_hitrost0; g=9.81, l=1.0)
    lastna_frekvenca = sqrt(g/l)
    return kot0 * cos(lastna_frekvenca * cas) +
           (kotna_hitrost0/lastna_frekvenca) * sin(lastna_frekvenca * cas)
end

# ------------------------------
# 2. Primerjava nelinearno vs linearno
# ------------------------------
zacetni_kot = 0.5
zacetna_kotna_hitrost = 0.0
cas, resitev = resi_nihalo(zacetni_kot, zacetna_kotna_hitrost, 10.0)

kot_linearno = [linearno_nihalo(ti, zacetni_kot, zacetna_kotna_hitrost) for ti in cas]

plot(cas, resitev[:,1], label="Matematično nihalo")
plot!(cas, kot_linearno, linestyle=:dash, label="Harmonično nihalo")
xlabel!("čas")
ylabel!("kot(t)")
savefig("primerjava_nihala.png")

# ------------------------------
# 3. Izračun periode matematičnega nihala
# ------------------------------
function izracunaj_periodo(zacetni_kot; zacetna_kotna_hitrost=0.0, g=9.81, l=1.0)
    cas, resitev = resi_nihalo(zacetni_kot, zacetna_kotna_hitrost, 30.0; epsilon=1e-8, sigma=0.9)
    spremembe = findall(i -> resitev[i,1] * resitev[i+1,1] < 0, 1:(length(cas)-1))
    if length(spremembe) >= 2
        return cas[spremembe[2]] - cas[spremembe[1]]
    else
        return NaN
    end
end

# ------------------------------
# 4. Graf nihajnega časa v odvisnosti od energije
# ------------------------------
zacetni_koti = [0.1, 0.5, 1.0, 1.5]  # radiani
periode = [izracunaj_periodo(k) for k in zacetni_koti]
energije = [energija_nihalo([k, 0.0]) for k in zacetni_koti]

plot(energije, periode, marker=:o, xlabel="Energija", ylabel="Perioda",
     label="Nihajni čas matematičnega nihala")
savefig("perioda_vs_energija.png")
