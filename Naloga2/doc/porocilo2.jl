#' # Ploščina zanke Bezierjeve krivulje
#' Tina Brdnik
#'
#' ## Opis naloge
#' Bézierjeve krivulje se uporabljajo za modeliranje gladkih in prilagodljivih krivulj.
#' Definirane so s pomočjo kontrolnih točk in Bernsteinovih polinomov. Krivulja stopnje $n$
#' z kontrolnimi točkami $P_0, P_1, \dots, P_n$ je podana parametrično kot:
#'
#' 
#' $B(t) = \sum_{i=0}^n \binom{n}{i}(1-t)^{n-i}t^i P_i, \quad t \in [0,1].$
#' 
#'
#' Ploščina zanke, ki jo omejuje Bézierjeva krivulja, se izračuna z integracijo diferencialne
#' ploščine:
#'
#' 
#' $A = \tfrac{1}{2}\int_0^1 \big(x(t)y'(t) - x'(t)y(t)\big)\,dt.$
#' 
#' ## Opis rešitve
#' Najprej smo implementirali funkcijo za računanje Bézierjeve krivulje na podlagi danih
#' kontrolnih točk, kjer smo za uteži uporabili Bernsteinove polinome in lastno implementacijo
#' binomskih koeficientov. Nato smo sestavili funkcijo za odvod krivulje, ki temelji na dejstvu,
#' da je odvod Bézierjeve krivulje stopnje $n$ spet Bézierjeva krivulja stopnje $n-1$ z novimi
#' kontrolnimi točkami $n(P_{i+1} - P_i)$.
#'
#' Integral $x(t)y'(t) - x'(t)y(t)$ smo nato integrirali z uporabo sestavljenega Simpsonovega
#' pravila. Na ta način smo dobili ploščino zanke, ki jo določa dana Bézierjeva krivulja.
#'
#' ## Rezultati
#' Kot primer smo obravnavali Bézierjevo krivuljo s kontrolnimi točkami
#'
#' $[(0,0), (1,1), (2,3), (1,4), (0,4), (-1,3), (0,1), (1,0)].$
#'
#' Na levi sliki je prikazana Bézierjeva krivulja (modra) skupaj s kontrolnimi točkami
#' (rdeče). Na desni sliki je prikazan integrand, ki ga uporabimo pri izračunu ploščine.
#' Dobljena ploščina zanke je enaka $1.9662004662$.
#'
include("demo.jl")
using .Naloga2
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
area, plot_bezier = demo_bezier(kontrolne_tocke; n=2000, nt=10_000)
display(plot_bezier)
