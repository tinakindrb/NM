# Matematično nihalo
### Avtor: Tina Brdnik

## Opis naloge
Naloga obravnava matematično nihalo, ki ga opišemo z nelinearno diferencialno enačbo drugega reda. Enačbo prevedemo na sistem prvega reda in rešujemo z numerično metodo Dormand–Prince 5. reda z adaptivnim korakom.  
Rešitve matematičnega nihala primerjamo s harmoničnim nihalom, kjer dobimo analitično rešitev. Poleg časovnega poteka odmika izračunamo tudi nihajni čas za različne začetne kote ter ga prikažemo glede na začetno energijo sistema.


## Primer uporabe
Primer uporabe je prikazan v programu, ki se nahaja v 'doc/demo.jl'. Poženemo ga v interaktivni zanki Julie z ukazi:

```jl
include("Naloga3/doc/demo.jl")
plt = demo_nihalo()
```

## Testi
Teste poženemo z ukazom:
```jl
include("Naloga3/test/runtests.jl")
```

## Poročilo PDF

Poročilo pripravimo z ukazoma:
```jl
using Weave

Weave.weave("Naloga3/doc/porocilo.jl",doctype="minted2pdf", out_path="Naloga3/pdf")
```
Poročilo se nahaja v Naloga3/pdf/porocilo.pdf

