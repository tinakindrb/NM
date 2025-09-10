#' # Inverzna potenčna metoda za zgornjo Hessenbergovo matriko
#' Tina Brdnik
#' 
#' ## Opis naloge
#'
#' V poročilu je opisana implementacija inverzne potenčne metode za določanje lastnih vektorjev matrike. 
#' Metoda temelji na zaporedju: 
#'
#' $x^{(n+1)} = \frac{A_{\lambda}^{-1} \, x^{(n)}}{\left\| A_{\lambda}^{-1} x^{(n)} \right\|_{\infty}}, \quad A_{\lambda} = A - \lambda I$
#' 
#' ki konvergira k lastnemu vektorju za lastno vrednost matrike
#' $A$ najbližji λ.
#' 
#' Za zmanjšanje števila potrebnih operacij, matriko $A$ najprej
#' s Householderjevimi zrcaljenji pretvorimo v zgornjo Hessenbergovo obliko. 
#' Razcep LU take matrike zahteva bistveno manj računskih operacij kot običajni LU razcep. Ker konjugacija s Hausholderjevimi zrcaljenji ohranja lastne vrednosti, ima tako dobljena zgornje Hessenberova matrika iste lastne vrednosti, kot originalna. 
#'
#' ## Opis rešitve
#'
#' Implementacija algoritma poteka v več korakih:
#'
#' 1. Pretvorba matrike v Hessenbergovo obliko  
#' Matriko $A$ najprej pretvorimo v zgornjo Hessenbergovo obliko s pomočjo
#' Householderjevih zrcaljenj. Na ta način odpravimo elemente pod drugo poddiagonalo.
#' Dobljeno matriko označimo z $H$. Zanjo velja $H = Q^T A Q$, kjer je $Q$
#' ortogonalna matrika. Ker imata $A$ in $H$ enake lastne vrednosti,
#' lahko nadaljnje izračune izvajamo kar na $H$.
#'
#' 2. LU razcep  
#' Za reševanje sistemov izvedemo LU-razcep zgornje Hessenbergove matrike.  
#' Prednost take oblike je v tem, da je razcep bistveno učinkovitejši.  
#' Zaradi redke strukture poteka eliminacija le ob diagonali, kar zagotavlja časovno zahtevnost $\mathcal{O}(n^2)$.  
#' Faktor $L$ je spodnje tridiagonalna matrika z enicami na diagonali, faktor $U$ pa zgornjetrikotna matrika.
#'
#' 3. Inverzna iteracija  
#' Inverzna potenčna metoda poišče lastno vrednost matrike, ki je najbližja
#' podanemu približku. Na vsakem koraku rešimo sistem $(A - \lambda I)y = x$
#' in dobljeni vektor normiramo. Ker račun poteka na Hessenbergovi obliki in
#' uporabljamo LU razcep, je metoda bistveno hitrejša kot v splošnem primeru.
#'
#' 4. Rekonstrukcija lastnega vektorja  
#' Če je $A = Q^T H Q$ in je $v$ lastni vektor matrike $H$, potem je $Qv$
#' lastni vektor matrike $A$ pri isti lastni vrednosti.
#'
#'
#' ## Primer uporabe
#'
#' Inverzno iteracijo smo preizkusili na polinomu
#'
#' $p(x) = x^3 - 2x^2 - 4x + 6.$
#'
#' Njegove ničle so enake lastnim vrednostim spremljevalne matrike, ki pripada
#' polinomu. Inverzno iteracijo smo izvedli na tej matriki.
#'
#' Na spodnji sliki je narisan graf polinoma skupaj z njegovimi ničlami.
#' Označeni so tudi začetni približki na realni osi in pripadajoče ničle,
#' proti katerim metoda konvergira. Tako dobimo območja konvergence
#' inverzne iteracije na realni osi.
include("demo.jl")
using .Naloga1
nariši_polinom_s_konvergenco([1.0, -2.0, -4.0, 6.0])








