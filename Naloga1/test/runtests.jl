using Naloga1, Test, LinearAlgebra

"""
TESTI funkcija hessenberg
"""

@testset "Hessenbergova redukcija" begin
    # Primer 1: naključna 4x4 matrika
    A = randn(4,4)
    H, Q = Naloga1.hessenberg(A)

    # (1) Hessenbergova oblika: elementi pod drugo poddiagonalo so ≈ 0
    for i = 3:4, j = 1:i-2
        @test isapprox(H.H[i,j], 0.0; atol=1e-12)
    end

    # (2) Podobnostna zveza
    @test Q' * A * Q ≈ H.H atol=1e-10

    # (3) Q ortogonalna
    @test Q * Q' ≈ I(4) atol=1e-12
    @test Q' * Q ≈ I(4) atol=1e-12


    # Primer 2: identiteta
    A = Matrix(I, 5, 5)
    H, Q = Naloga1.hessenberg(A)
    @test H.H ≈ A atol=1e-12
    @test Q * Q' ≈ I(5) atol=1e-12


    # Primer 3: že Hessenbergova matrika
    A = [1.0 2.0 3.0;
         4.0 5.0 6.0;
         0.0 7.0 8.0]
    H, Q = Naloga1.hessenberg(A)
    @test H.H ≈ A atol=1e-12
    @test Q * Q' ≈ I(3) atol=1e-12
end

"""
TESTI podatkovni tip SpTridiag
"""

@testset "Podatkovni tip SpTridiag je definiran" begin
    d = [1.0, 1.0, 1.0, 1.0]
    sd = [2.0, 3.0, 4.0]
    L = SpTridiag(d, sd)
    @test L.d == d
    @test L.sd == sd
end

@testset "Rešitev sistema Lx=b z L\\b" begin
    d = [1.0, 1.0, 1.0, 1.0]
    sd = [2.0, 3.0, 4.0]
    L = SpTridiag(d, sd)

    # Sestavimo znano rešitev
    x_exact = [1.0, 2.0, 3.0, 4.0]
    b = zeros(4)
    b[1] = d[1]*x_exact[1]
    for i = 2:4
        b[i] = sd[i-1]*x_exact[i-1] + d[i]*x_exact[i]
    end

    # Preverimo reševanje
    x_computed = L \ b
    @test x_computed ≈ x_exact
end

@testset "Napaka ob napačni dimenziji" begin
    d = [1.0, 1.0]
    sd = [2.0]
    L = SpTridiag(d, sd)
    b = [1.0, 2.0, 3.0]   # predolg vektor
    @test_throws DimensionMismatch L \ b
end

"""
TESTI podatkovni tip ZgornjiHessenberg
"""

# pomožna funkcija za razširitev SpTridiag v polno matriko (za testiranje)
function Base.Matrix(L::SpTridiag)
    n = length(L.d)
    M = Matrix{eltype(L.d)}(I, n, n)   # diagonala = enice
    for i = 1:n-1
        M[i+1, i] = L.sd[i]            # poddiagonala iz L
    end
    return M
end

@testset "Konstruktor ZgornjiHessenberg" begin
    H = [-149.0   42.2037  -156.3165;
         -537.6783 152.5511 -554.9272;
           0.0      0.0728    2.4489]
    zh = ZgornjiHessenberg(H)
    @test zh.H ≈ H atol=1e-8

    I = Matrix(LinearAlgebra.I, 3, 3)
    zh = ZgornjiHessenberg(I)
    @test zh.H ≈ I atol=1e-12

    H2 = [1.0    -3.5970   -0.2481;
          -8.0623  14.0462   2.8308;
           0.0     0.8308  -0.0462]
    zh = ZgornjiHessenberg(H2)
    @test zh.H ≈ H2 atol=1e-8
end

@testset "LU razcep zgornje Hessenbergove matrike" begin
    # Primer 1
    H = [-149.0 -50.0 -154.0;
          537.0 180.0  546.0;
          -27.0  -9.0  -25.0]
    zh = ZgornjiHessenberg(Naloga1.hessenberg(H))
    L, U = lu(zh)
    @test Matrix(L)*U ≈ zh.H atol=1e-10

    # Primer 2
    A = [4.0 3.0 2.0;
         2.0 3.0 1.0;
         1.0 1.0 4.0]
    zh = ZgornjiHessenberg(A)
    L, U = lu(zh)
    @test Matrix(L)*U ≈ zh.H atol=1e-10

    # Primer 3
    B = [2.0 -1.0  0.0;
         -1.0  2.0 -1.0;
          0.0 -1.0  2.0]
    zh = ZgornjiHessenberg(B)
    L, U = lu(zh)
    @test Matrix(L)*U ≈ zh.H atol=1e-10
end



