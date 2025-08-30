# Inverzna potenčna metoda za zgornjo Hessenbergovo matriko
### Avtor: Tina Brdnik

## Opis naloge
Naloga implementira inverzno potenčno metodo za izračun lastnih vrednosti in lastnih vektorjev matrike $A$.  
Za pospešitev izračuna se matrika najprej pretvori v zgornjo Hessenbergovo obliko s Householderjevimi zrcaljenji, nato pa se uporabi optimiziran LU razcep.  
Metoda je preizkušena tudi na matriki polinoma, kjer lastne vrednosti ustrezajo ničlam polinoma.

## Primer uporabe
Primer uporabe je prikazan v programu, ki se nahaja v 'doc/demo.jl'. Poženemo ga v interaktivni zanki Julie z ukazi:

```jl
include("Naloga1/doc/demo.jl")
using Naloga1
nariši_polinom_s_konvergenco([koeficienti polinoma])
```

## Testi
Teste poženemo z ukazom:
```jl
include("Naloga1/test/runtests.jl")
```

## Poročilo PDF

Poročilo pripravimo z ukazoma:
```jl
using Weave

Weave.weave("Vaja01/doc/01uvod.jl",doctype="minted2pdf", out_path="Vaja01/pdf")
```
Poročilo se nahaja v Naloga1/pdf/porocilo.pdf

