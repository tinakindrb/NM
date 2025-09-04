using Test
include("../src/Naloga2.jl")
using .Naloga2

"""
TESTI za Simpsonovo pravilo
"""

@testset "Simpsonovo pravilo" begin
    # Integral x^2 od 0 do 1 = 1/3
    f1(x) = x^2
    result1 = Naloga2.simpsonovo_pravilo(f1, 0.0, 1.0, 4)
    @test isapprox(result1, 1/3, atol=1e-10)

    # Integral sin(x) od 0 do π = 2
    f2(x) = sin(x)
    result2 = Naloga2.simpsonovo_pravilo(f2, 0.0, Float64(π), 1000)
    @test isapprox(result2, 2.0, atol=1e-10)

    # Test napake, če je n liho
    @test_throws ErrorException Naloga2.simpsonovo_pravilo(f1, 0.0, 1.0, 3)
end

"""
TESTI za porazdelitveno funkcijo standardne normalne porazdelitve
"""

@testset "Standardna normalna CDF F(x)" begin
    # F(0) = 0.5
    @test isapprox(F(0.0), 0.5, atol=1e-12)

    # Simetrija: F(-x) ≈ 1 - F(x)
    @test isapprox(F(-1.5), 1 - F(1.5), atol=1e-12)

    # F(1) ≈ 0.841344746...
    @test isapprox(F(1.0), 0.841344746, atol=1e-10)

    # Negativen x, recimo -2
    @test isapprox(F(-2.0), 0.0228, atol=1e-3)

    # Zelo majhen x (rep na levi strani)
    @test Naloga2.F(-10.0) ≈ 0.0 atol=1e-15

    # Zelo velik x (rep na desni strani)
    @test Naloga2.F(10.0) ≈ 1.0 atol=1e-10
end


@testset "Binomski koeficient" begin
    # Znane vrednosti
    @test Naloga2.binomski(5,0) == 1
    @test Naloga2.binomski(5,1) == 5
    @test Naloga2.binomski(5,2) == 10
    @test Naloga2.binomski(5,5) == 1
    # Vrednosti zunaj obsega
    @test Naloga2.binomski(5,-1) == 0
    @test Naloga2.binomski(5,6) == 0
end


@testset "Bézierjeva krivulja" begin
    # Preprost primer: kontrolne točke (0,0) in (1,1) -> premica y=x
    kontrolne_tocke = [0.0 0.0; 1.0 1.0]
    t = [0.0, 0.5, 1.0]
    rezultat = Naloga2.bezier(kontrolne_tocke, t)
    @test isapprox(rezultat, [0.0 0.0; 0.5 0.5; 1.0 1.0], atol=1e-12)

    # Lastnost: na začetku in koncu krivulja gre skozi prvi in zadnji kontrolni točki
    kontrolne_tocke2 = [0.0 0.0; 1.0 2.0; 2.0 0.0]
    krivulja = Naloga2.bezier(kontrolne_tocke2, [0.0, 1.0])
    @test krivulja[1,:] == kontrolne_tocke2[1,:]
    @test krivulja[end,:] == kontrolne_tocke2[end,:]
end

@testset "Bézierjeva krivulja z odvodom" begin
    # linearna Bézierjeva krivulja med (0,0) in (1,1) -> premica y=x
    kontrolne_tocke = [0.0 0.0; 1.0 1.0]
    t = [0.0, 0.5, 1.0]
    krivulja, odvodi = Naloga2.bezier_in_odvod(kontrolne_tocke, t)

    # Krivulja mora biti na premici y=x
    @test all(isapprox.(krivulja[:,1], krivulja[:,2]; atol=1e-12))

    # Odvod mora biti konstanten (1,1)
    for i in 1:length(t)
        @test isapprox(odvodi[i,1], 1.0; atol=1e-12)
        @test isapprox(odvodi[i,2], 1.0; atol=1e-12)
    end

    # kvadratna Bézierjeva krivulja s simetrijo
    kontrolne_tocke2 = [0.0 0.0; 1.0 2.0; 2.0 0.0]
    t2 = [0.0, 0.5, 1.0]
    krivulja2, odvodi2 = Naloga2.bezier_in_odvod(kontrolne_tocke2, t2)

    # Na začetku in koncu se odvod ujema s tangentami na poligon
    # Odvod na začetku: n*(P1-P0) = 2*(1,2)-(0,0) = (2,4)
    @test isapprox(odvodi2[1,1], 2.0; atol=1e-12)
    @test isapprox(odvodi2[1,2], 4.0; atol=1e-12)

    # Odvod na koncu: n*(Pn-P{n-1}) = 2*((2,0)-(1,2)) = (2,-4)
    @test isapprox(odvodi2[end,1], 2.0; atol=1e-12)
    @test isapprox(odvodi2[end,2], -4.0; atol=1e-12)
end

@testset "Ploščina Bézierjeve krivulje" begin
    # Lastnost: ploščina mora biti končno število
    kontrolne_tocke = [0.0 0.0; 1.0 1.0; 2.0 0.0]
    A = Naloga2.ploscina_bezier(kontrolne_tocke; n=1000)
    @test isfinite(A)

    # Lastnost: če krivuljo obrneš, se ploščina spremeni za predznak
    kontrolne_tocke_obrnjene = reverse(kontrolne_tocke, dims=1)
    A2 = Naloga2.ploscina_bezier(kontrolne_tocke_obrnjene; n=1000)
    @test isapprox(abs(A), abs(A2); atol=1e-6)

    # Lastnost: ploščina mora biti blizu ničle za "degenerirane" primere
    kontrolne_tocke_premica = [0.0 0.0; 0.5 0.5; 1.0 1.0]
    A3 = Naloga2.ploscina_bezier(kontrolne_tocke_premica; n=1000)
    @test isapprox(A3, 0.0; atol=1e-6)

    kontrolne_tocke2 = [
        0.0 0.0;
        1.0 1.0;
        2.0 3.0;
        1.0 4.0;
        0.0 4.0;
       -1.0 3.0;
        0.0 1.0;
        1.0 0.0
    ]
    A2 = Naloga2.ploscina_bezier(kontrolne_tocke2; n=2000)
    @test isapprox(A2, 1.966200466200; atol=1e-10)
end

