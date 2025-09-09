#' # Matematično nihalo
#' Tina Brdnik
#'
#' ## Opis naloge
#' Naloga obravnava matematično nihalo, katerega odmik $kot(t)$ je podan z
#' diferencialno enačbo:
#' $$
#' \ddot{\theta}(t) + \frac{g}{l}\sin(\theta(t)) = 0, \quad \theta(0) = \theta_0, \quad \dot{\theta}(0) = \omega_0.
#' $$
#' Enačbo drugega reda prevedemo na sistem prvega reda:
#' $$
#' \begin{cases}
#' \theta'(t) = \omega(t), \\
#' \omega'(t) = -\dfrac{g}{l}\sin(\theta(t)).
#' \end{cases}
#' $$
#' Za različne začetne pogoje primerjamo rešitev z nihanjem harmoničnega nihala,
#' ki je podano z linearizirano enačbo:
#' $$
#' \ddot{\theta}(t) + \frac{g}{l}\theta(t) = 0.
#' $$
#' Pri harmoničnem nihalu je nihajni čas neodvisen od začetnega odmika, medtem
#' ko pri matematičnem nihalu narašča z večjo začetno energijo.
#' ## Opis rešitve
#'
#' Implementacija je razdeljena na več funkcij, združenih v modulu `Naloga3`:
#'
#' 1. **Metoda DOPRI5**  
#' Funkcija `dopri5` je lastna implementacija metode Dormand–Prince 5. reda (DOPRI5) z
#' adaptivnim korakom. Metoda uporablja oceno napake med približkom 5. in 4. reda ter
#' na tej osnovi povečuje ali zmanjšuje dolžino koraka. S tem zagotavlja visoko
#' natančnost pri optimalnem številu korakov.
#'
#' 2. **Sistem matematičnega nihala**  
#' Enačbo drugega reda smo pretvorili v sistem prvega reda in zapisali funkcijo
#' `sistem_nihalo(t,y)`, ki vrača odvode $[\theta', \omega']$ za dano stanje.
#'
#' 3. **Energija sistema**  
#' Funkcija `energija_nihalo(y)` izračuna mehansko energijo nihala kot:
#'
#' $$
#' E = \tfrac{1}{2}\,\omega^2 + \frac{g}{l}\,(1 - \cos(\theta)).
#' $$
#'
#' Energija se v točni rešitvi ohranja, zato je uporabna za preverjanje pravilnosti
#' numerične metode.
#'
#' 4. **Reševanje nihala**  
#' Funkcija `resi_nihalo(zacetni_kot, zacetna_kotna_hitrost, cas_konec)` uporabi
#' metodo `dopri5` in funkcijo `sistem_nihalo` za izračun rešitve na podanem časovnem
#' intervalu. Rezultat sta vektor časov in matrika, ki vsebuje vrednosti kota in
#' kotne hitrosti v teh časih.
