#' # Porazdelitvena funkcija za N(0, 1)
#' Tina Brdnik
#'
#' ## Opis naloge
#' Naloga obravnava porazdelitveno funkcijo standardne normalne porazdelitve N(0,1), ki je definirana kot:
#'
#' ```math
#' F(x) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^x e^{-t^2/2}\,dt
#' ```
#'
#' Ta funkcija podaja verjetnost, da je slučajna spremenljivka z normalno porazdelitvijo manjša ali enaka vrednosti x.
#' ## Opis rešitve
#'
#' Implementacija algoritma poteka v več korakih:
#'
#' 1. Numerična integracija z Simpsonovim pravilom  
#' Uporabimo sestavljeno Simpsonovo pravilo za izračun določenega integrala
#' funkcije gostote normalne porazdelitve. To pravilo združuje vrednosti funkcije
#' v krajiščih in središčih podintervalov ter zagotavlja visoko natančnost.
#'
#' 2. Posebni primeri  
#' Za $x = 0$ lahko rezultat podamo neposredno kot $F(0) = 0.5$.
#' Za zelo velike vrednosti $|x|$ uporabimo asimptotsko aproksimacijo repa
#' normalne porazdelitve, da se izognemo nepotrebnemu računanju na daljših
#' intervalih.
#'
#' 3. Simetrija funkcije  
#' Ker velja $F(-x) = 1 - F(x)$, za negativne vrednosti $x$ izračunamo
#' rezultat prek pozitivnega argumenta. Na ta način izkoristimo lastnosti
#' normalne porazdelitve in poenostavimo izračun.
#'
#' 4. Kombinacija metod  
#' Za $x > 0$ integral izračunamo z uporabo Simpsonovega pravila na intervalu
#' $[0,x]$, pri tem pa dodamo začetno vrednost $0.5$, ki ustreza $F(0)$.
#' Tako dobimo splošno funkcijo $F(x)$, ki vrača približek porazdelitvene
#' funkcije standardne normale.

#' ## Rezultati
#'
#' Funkcijo $F(x)$ smo primerjali z referenčno rešitvijo, izračunano s pomočjo
#' funkcije 'erf'. Na spodnjem grafu sta prikazani obe krivulji, ki se praktično
#' popolnoma prekrivata. Izračunali smo tudi delež točk, kjer je bila napaka manjša
#' od dovoljene tolerance, in dobili rezultat 100 %.
include("demo.jl")
using .Naloga2
plt, pct = demo_F()
display(plt)