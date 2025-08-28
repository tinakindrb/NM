module Naloga1

using LinearAlgebra

export ZgornjiHessenberg, hessenberg

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


end # modul Naloga1
