# Porazdelitvena funkcija za N(0, 1) in Ploščina zanke Bezierjeve krivulje
### Avtor: Tina Brdnik

## Opis naloge
Naloga implementira dva numerična algoritma:  

1. **Porazdelitvena funkcija standardne normalne porazdelitve**  
Funkcija $F(x)$ je definirana kot integral gostote normalne porazdelitve.  
Za izračun smo uporabili sestavljeno Simpsonovo pravilo, pri čemer smo pri
negativnih vrednostih izkoristili simetrijo $F(-x) = 1 - F(x)$, za zelo
velike vrednosti pa asimptotsko aproksimacijo repa porazdelitve.  

2. **Ploščina zanke Bézierjeve krivulje**  
Bézierjeva krivulja je definirana s kontrolnimi točkami in Bernsteinovimi
polinomi. Odvod krivulje smo izračunali analitično, saj je odvod Bézierjeve
krivulje spet Bézierjeva krivulja nižje stopnje. Ploščina je nato
izračunana po formuli
$A = \tfrac{1}{2}\int_0^1 \big(x(t)y'(t) - x'(t)y(t)\big)\,dt$

pri čemer smo integracijo izvedli z uporabo Simpsonovega pravila.  


## Primer uporabe
Primera uporabe sta prikazana v programu, ki se nahaja v 'doc/demo.jl'. Poženemo ga v interaktivni zanki Julie z ukazi:

**Za prvo nalogo**
```jl
include("Naloga2/doc/demo.jl")
plt, pct = demo_F()
display(plt)
```
**Za drugo nalogo**
```jl
include("Naloga2/doc/demo.jl")
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
```

## Testi
Teste poženemo z ukazom:
```jl
include("Naloga2/test/runtests.jl")
```

## Poročilo PDF

Poročili pripravimo z ukazia:
```jl
using Weave

Weave.weave("Naloga2/doc/porocilo1.jl",doctype="minted2pdf", out_path="Naloga2/pdf")

Weave.weave("Naloga2/doc/porocilo2.jl",doctype="minted2pdf", out_path="Naloga2/pdf")
```
Poročili se nahajata v Naloga2/pdf pod porocilo1.pdf in porocilo2.pdf

