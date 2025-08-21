module Naloga1

using LinearAlgebra

export SpTridiag, ZgornjiHessenberg

struct SpTridiag
    d   # glavna diagonala (dolžine n)
    sd  # poddiagonala (dolžine n-1)
end

import Base: \

"""
L \\ b
Reši sistem Lx = b, kjer je L spodnja tridiagonalna matrika tipa SpTridiag.
"""
function \(L::SpTridiag, b::Vector)
    n = length(L.d)
    if length(b) != n
        throw(DimensionMismatch("Dimenzija vektorja b se ne ujema z matriko L"))
    end
    x = similar(b)
    x[1] = b[1] / L.d[1]
    for i = 2:n
        x[i] = (b[i] - L.sd[i-1]*x[i-1]) / L.d[i]
    end
    return x
end

struct ZgornjiHessenberg
    H   # kvadratna zgornja Hessenbergova matrika
end

"""
lu(Hh::ZgornjiHessenberg)

Izvede LU razcep zgornje Hessenbergove matrike H brez pivotiranja.
Vrne (L::SpTridiag, U::Matrix).
"""
function lu(Hh::ZgornjiHessenberg)
    H = Hh.H
    n = size(H, 1)

    Lsd = zeros(eltype(H), n-1)
    U = zeros(eltype(H), n, n)

    # prva vrstica U
    U[1, :] .= H[1, :]

    for i = 1:n-1
        # multiplikator
        Lsd[i] = H[i+1, i] / U[i, i]

        # vrstica i+1 v U od stolpca i+1 dalje
        for j = i+1:n
            U[i+1, j] = H[i+1, j] - Lsd[i] * U[i, j]
        end
    end

    L = SpTridiag(ones(eltype(H), n), Lsd)
    return L, U
end

function hessenberg(A::Matrix)
    n = size(A, 1)
    Q = Matrix{eltype(A)}(I, n, n)
    H = copy(A)

    for k = 1:n-2
        x = H[k+1:n, k]

        # če je stolpec že v Hessenberg obliki, preskoči
        if norm(x[2:end]) < 1e-14
            continue
        end

        e = zeros(eltype(A), length(x))
        e[1] = norm(x)

        if norm(x) == 0
            u = zeros(eltype(A), length(x))
        else
            u = sign(x[1]) * e + x
            u /= norm(u)
        end

        P = Matrix{eltype(A)}(I, n-k, n-k) - 2 * (u * u')

        H[k+1:n, k:n] = P * H[k+1:n, k:n]
        H[1:n, k+1:n] = H[1:n, k+1:n] * P
        Q[k+1:n, :]   = P * Q[k+1:n, :]
    end

    return ZgornjiHessenberg(H), Q'
end


end # modul Naloga1
