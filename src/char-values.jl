export charλ,
       charA,
    charB,
    ce,
    se
"""
charλ(q,ν,k)

char value λ_(ν+k) for Mathieu's equation

y'' + (λ_(ν+k) - 2 q cos( 2z )) y = 0

where

q ∈ ℝ       - parameter
ν ∈ [-1,1]  - fractional part of the non-integer order
k ∈ ℤ⁺      - range of integer parts of the order

"""
function charλ(q::Real, nu_::Real; k::UnitRange=1:1) # reduced = true
    #nu = reduced ? rem(nu_+1,2)-1 : nu_;
    nu = rem(nu_+1,2)-1;

    # Set matrix size using formula from Shirts paper (1993), Eqs. (2.1)-(2.2).
    nu0 = nu + maximum(k)
    C = (8.46 + 0.444*nu0)/(1 + 0.085*nu0)
    D =  (0.24 + 0.0214*nu0)/(1 + 0.059*nu0)
    N = ceil(Int, (nu0 + 2 + C*abs(q)^D)/2) # matrix size is 2N+1

    (two, q2, nu2) = float.(promote(2, q, nu)) # d0 and d1 must be of the same type
    d0 = (two .* (-N:N) .- nu2).^2
    d1 = fill(q2, 2 * N)
    A = SymTridiagonal(d0, d1)
    a = eigvals!(A, k)
    return a
end

"""
charA(q; k=0:4)

char value A_k for Mathieu's equation

y'' + (A_k - 2 q cos( 2z )) y = 0

where

q ∈ ℝ  - parameter
k ∈ ℤ⁺ - eigenvalue index

"""
function charA(q::Real; k::UnitRange=1:1)
    all(x -> x >= 0, k) || throw(DomainError(k, "Indices must be non-negative integers."))

    # Boolean indices of even and odd n values
    ie = map(iseven, k)
    io = map(!, ie)

    a = Array{Float64}(undef ,length(k))
    k1 = k .+ 1
    a[ie] = charλ(abs(q), 0.0; k = k1)[ie]
    if q>=0
        a[io] = charλ(q, one(q); k = k1)[io]
    else
        if 0 in k # maybe not the cleanest way to do it
            a[io] = charλ(abs(q), one(q); k = k[2]:last(k))[io[2:end]]
        else
            a[io] = charλ(abs(q), one(q); k=k)[io]
        end
    end
    return a
end


"""
charB(q,k)

char value B_k for Mathieu's equation

y'' + (B_k - 2 q cos( 2z )) y = 0

where

q ∈ ℝ  - parameter
k ∈ ℤ  - eigenvalueindex

"""
function charB(q::Real; k::UnitRange=1:1)
    all(x -> x > 0, k) || throw(DomainError(k, "Indices must be positive integers."))
    # Boolean indices of even and odd n values
    ie = map(iseven, k)
    io = map(!, ie)

    b = Array{Float64}(undef, length(k))
    b[ie] = charλ(q,0.0,k=k)[ie]
    if q>=0
        b[io] = charλ(q,1.0,k=k)[io]
    else
        b[io] = charλ(abs(q),1.0,k = (k .+ 1))[io]
    end
    return b
end

#Mathieu cosine function:
function ce(m::Int64, q::Real, x::AbstractVector{Float64})
    accuracy = 0.00001
    if m%2 == 0
        a = charA(q, k=m:m)[1]
        A = [0., 1., 0., a/q, 0., ((a-4.)*a/(q*q))-2., 0.]
        i = 8
        #Normalization as each new coefficient is computed:
        N = 2*A[2]^2 + sum(A[4:end].^2)
        A = A ./sqrt(N)
        while (abs(A[i-2])>accuracy)
            push!(A, ((a-(i-4.)*(i-4.))/q)*A[i-2]- A[i-4])
            push!(A, 0.)
            N = 2*A[2]^2 + sum(A[4:end].^2)
            A = A ./sqrt(N)
            i += 2
        end
        t = zeros(length(x))
        for j in 2:length(A)
            t += A[j]*cos.((j-2)*x)
        end
        return t
    end
    if m%2 == 1
        a = charA(q, k=m:m)[1]
        A = [1., 0., (a-q-1.)/q, 0., ((a-9.)/q)*((a-q-1.)/q)-1., 0., ((a-25.)/q)*(((a-9.)/q)*((a-q-1.)/q)-1.)-(a-q-1.)/q, 0.]
        i = 9
        #Normalization as each new coefficient is computed:
        N = sum(A[1:end].^2)
        A = A ./sqrt(N)
        while (abs(A[i-2])>accuracy)
            push!(A, ((a-(i-2)*(i-2))/q)*A[i-2] - A[i-4])
            push!(A, 0.)
            N = sum(A[1:end].^2)
            A = A ./sqrt(N)
            i += 2
        end
        t = zeros(length(x))
        for j in 1:length(A)
            t += A[j]*cos.((j)*x)
        end
        return t
    end
end

#Mathieu sine function:
function se(m::Int64, q::Real, x::AbstractVector{Float64})
    accuracy = 0.00001
    if m%2 == 0
        a = charB(q, k=m:m)[1]
        B = [0., 1., 0., (a-4.)/q, 0.]
        i = 6
        #Normalization as each new coefficient is computed:
        N = sum(B[2:end].^2)
        B = B ./sqrt(N)
        while (abs(B[i-2])>accuracy)
            push!(B, ((a-(i-2.)*(i-2.))/q)*B[i-2]- B[i-4])
            push!(B, 0.)
            N = sum(B[2:end].^2)
            B = B ./sqrt(N)
            i += 2
        end
        t = zeros(length(x))
        for j in 2:length(B)
            t += B[j]*sin.((j)*x)
        end
        return t
    end
    if m%2 == 1
        a = charB(q, k=m:m)[1]
        B = [1., 0., (a+q-1.)/q, 0.]
        i = 5
        #Normalization as each new coefficient is computed:
        N = sum(B[1:end].^2)
        B = B ./sqrt(N)
        while (abs(B[i-2])>accuracy)
            push!(B, ((a-(i-2)*(i-2))/q)*B[i-2] - B[i-4])
            push!(B, 0.)
            N = sum(B[1:end].^2)
            B = B ./sqrt(N)
            i += 2
        end
        t = zeros(length(x))
        for j in 1:length(B)
            t += B[j]*sin.((j)*x)
        end
        return t
    end
end
