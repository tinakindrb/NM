module Naloga1

using LinearAlgebra

export ZgornjiHessenberg, hessenberg, SpTridiag, lu

"""
ZgornjiHessenberg(H)

Podatkovni tip za zgornje Hessenbergove matrike. Shranjuje kvadratno matriko H,
ki ima lahko neničelne elemente na glavni diagonali, na prvi poddiagonali in
na vseh naddiagonalnih mestih. Vsi elementi pod prvo poddiagonalo morajo biti 0.
"""
struct ZgornjiHessenberg
    H::Matrix{Float64}

    function ZgornjiHessenberg(H::AbstractMatrix)
        n, m = size(H)
        if n != m
            error("Matrika mora biti kvadratna.")
        end
        # preverimo Hessenberg lastnost: ničle pod drugo poddiagonalo
        for i = 3:n
            for j = 1:(i-2)
                if abs(H[i, j]) > 1e-12
                    error("Matrika ni zgornja Hessenbergova (element H[$i,$j] ≠ 0).")
                end
            end
        end
        new(Matrix{Float64}(H))
    end
end

Base.show(io::IO, Z::ZgornjiHessenberg) = print(io, "ZgornjiHessenberg:\n", Z.H)

"""
hessenberg(A)

Pretvori kvadratno matriko A v zgornjo Hessenbergovo obliko s pomočjo
Householderjevih zrcaljenj.

Rezultat je matrika H tipa ZgornjiHessenberg in ortogonalna matrika Q,
tako da velja: H = Q'*A*Q

# Vhod:
- A : kvadratna matrika

# Izhod:
- H : ZgornjiHessenberg (zgornja Hessenbergova matrika)
- Q : ortogonalna matrika
"""

function hessenberg(A::AbstractMatrix)
    n, m = size(A)

    H = Matrix{Float64}(A) # kopija A
    Q = Matrix{Float64}(I, n, n)

    for k = 1:n-2
        # vektor pod diagonalo v stolpcu k
        x = H[k+1:n, k]

        # Householderjev vektor
        e = zeros(length(x))
        e[1] = norm(x)
        if norm(x) == 0
            v = zeros(size(x))
        else
            v = sign(x[1]) * e + x
            v /= norm(v)
        end

        # refleksija
        P = Matrix{Float64}(I, length(v), length(v)) - 2 * (v * v')

        # posodobimo H in Q (uporabimo P na pravem podbloku)
        H[k+1:n, k:n] = P * H[k+1:n, k:n]
        H[1:n, k+1:n] = H[1:n, k+1:n] * P
        Q[k+1:n, :]   = P * Q[k+1:n, :]
    end

    #print(H)

    return ZgornjiHessenberg(H), Q'
end

"""
SpTridiag(sd)

Podatkovni tip za spodnje tridiagonalne matrike z enicami na glavni diagonali.

- sd : vektor dolžine n-1 z elementi poddiagonale
- glavna diagonala je privzeto [1, 1, ..., 1]
- vsi ostali elementi so 0
"""
struct SpTridiag
    d::Vector{Float64}
    sd::Vector{Float64}
end

"""
Matrix(L)

Pretvori elemente iz SpTridiag v polno matriko

- Glavna diagonala je podana z vektorjem d, 
- Poddiagonala je podana z vektorjem sd,
- Vsi ostali elementi matrike so enaki 0.

"""

function Matrix(L::SpTridiag)
    n = length(L.d)
    M = zeros(Float64, n, n)
    for i = 1:n
        M[i,i] = L.d[i] #glavna diagonala
    end
    for i = 1:n-1
        M[i+1,i] = L.sd[i] #poddiagonala
    end
    return M
end

"""
\\(L, b)

Reši sistem linearnih enačb L * x = b, kjer je L spodnja tridiagonalna
matrika tipa SpTridiag.

- Argument b je desna stran sistema.

Vrne rešitev x kot vektor dolžine n.
"""

import Base: \

function \(L::SpTridiag, b::Vector{Float64})
    n = length(L.d)
    if length(b) != n
        throw(ArgumentError("Dimenzije se ne ujemajo"))
    end

    x = similar(b)
    x[1] = b[1] / L.d[1]
    for i = 2:n
        x[i] = (b[i] - L.sd[i-1] * x[i-1]) / L.d[i]
    end
    return x
end

"""
lu(H)

Izvede LU razcep zgornje Hessenbergove matrike H.

Vrne:
- L : spodnja tridiagonalna matrika tipa SpTridiag z enicami na diagonali,
- U : zgornje trikotna matrika.
"""
function lu(Z::ZgornjiHessenberg)
    H = copy(Z.H)
    n = size(H, 1)

    d  = ones(Float64, n) #glavna diagonala L (same enice)
    sd = zeros(Float64, n-1) #poddiagonala L
    U  = zeros(Float64, n, n) #U kot polna matrika

    #inicializacija prve vrstice U
    U[1, :] = H[1, :]

    for i = 1:n-1
        #element L[i+1, i]
        sd[i] = H[i+1, i] / U[i, i]

        #naslednja vrstica v U
        U[i+1, i+1:n] = H[i+1, i+1:n] - sd[i] * U[i, i+1:n]
    end

    L = SpTridiag(d, sd)
    return L, U
end


end #modul Naloga1
