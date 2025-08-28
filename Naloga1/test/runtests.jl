using Naloga1, Test, LinearAlgebra

"""
TESTI podatkovni tip ZgornjiHessenberg
"""

@testset "ZgornjiHessenberg konstruktor" begin
    # Primer 1: konkretna Hessenbergova matrika
    H1 = [2.5   -1.2   3.0;
          -4.7   5.3  -6.1;
           0.0   7.4   8.8]
    zh1 = ZgornjiHessenberg(H1)
    @test zh1.H ≈ H1 atol=1e-8

    # Primer 2: identiteta
    I4 = Matrix{Float64}(I, 4, 4)
    zh2 = ZgornjiHessenberg(I4)
    @test zh2.H ≈ I4 atol=1e-12

    # Primer 3: naključna Hessenbergova (0 pod drugo poddiagonalo)
    H3 = [1.0    -2.0   3.5   0.0;
          5.1     0.7  -1.4   2.2;
          0.0    -3.3   4.0  -0.8;
          0.0     0.0   2.1   1.5]
    zh3 = ZgornjiHessenberg(H3)
    @test zh3.H ≈ H3 atol=1e-8

    # Primer 4: neveljavna matrika (mora sprožiti error)
    Hbad = [1.0  2.0  3.0;
            4.0  5.0  6.0;
            9.0  7.0  8.0]   # element (3,1) ≠ 0
    @test_throws ErrorException ZgornjiHessenberg(Hbad)

    # Primer 5: ne-kvadratna matrika (mora sprožiti error)
    Hrect = [1.0 2.0 3.0;
             4.0 5.0 6.0]   # 2x3 matrika
    @test_throws ErrorException ZgornjiHessenberg(Hrect)
end

"""
TESTI funkcija hessenberg
"""

@testset "hessenberg redukcija" begin
    # Primer 1: konkretna matrika
    A1 = [-149.0  -50.0  -154.0;
           537.0   180.0   546.0;
           -27.0    -9.0   -25.0]

    H1, Q1 = Naloga1.hessenberg(A1)

    # preveri rekonstrukcijo
    @test A1 ≈ Q1 * H1.H * Q1' atol=1e-10

    # Primer 2: identiteta
    A2 = Matrix{Float64}(I, 3, 3)
    H2, Q2 = Naloga1.hessenberg(A2)

    @test A2 ≈ Q2 * H2.H * Q2' atol=1e-10

    # Primer 3: "naključna" matrika
    A3 = [1.0  2.0  3.0;
          4.0  5.0  6.0;
          7.0  8.0  9.0]

    H3, Q3 = Naloga1.hessenberg(A3)

    @test A3 ≈ Q3 * H3.H * Q3' atol=1e-10
end



