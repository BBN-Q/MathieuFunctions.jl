export charλ,
       charA,
       charB
"""
charλ(q,ν,k)

char value λ_(ν+k) for Mathieu's equation

y'' + (λ_(ν+k) - 2 q cos( 2z )) y = 0

where

q ∈ ℝ       - parameter
ν ∈ [-1,1]  - fractional part of the non-integer order
k ∈ ℤ⁺      - range of integer parts of the order

"""

function charλ(q::Real, k_::StepRangeLen=1.0:1.0) # reduced = true
    @assert isinteger(step(k_)) && round(Int,step(k_))==1 "Range of eigenvalue indices must have unit step"

    k = floor(Int,first(k_)):floor(Int,last(k_));
    nu_ = first(k_)-floor(Int,last(k_));
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

function charλ(q::Real, nu::Number)
    charλ(q, nu:nu)
end

"""
char value A_k for Mathieu's equation

y'' + (A_k - 2 q cos( 2z )) y = 0

where

q ∈ ℝ  - parameter
k ∈ ℤ⁺ - eigenvalue index, or range of indices

"""

function charA(q::Real, k::StepRangeLen{<:Real}=1.0:1.0)
    all(x -> x >= 0, k) || throw(DomainError(k, "Indices must be non-negative integers."))

    # Boolean indices of even and odd n values
    k_ = map(x -> round(Int, x), k);
    ie = map(iseven, k_);
    io = map(!, ie);

    a = Array{Float64}(undef ,length(k));
    k1 = k .+ 1;
    a[ie] = charλ(abs(q), k1)[ie];
    if q>=0
        a[io] = charλ(q, k1)[io];
    else
        if 0 in k # maybe not the cleanest way to do it
            a[io] = charλ(abs(q), k[2]:last(k))[io[2:end]];
        else
            a[io] = charλ(abs(q), k)[io];
        end
    end
    return a
end

function charA(q::Real, k::Number)
    charA(q,k:k)
end

"""
charB(q,k)

char value B_k for Mathieu's equation

y'' + (B_k - 2 q cos( 2z )) y = 0

where

q ∈ ℝ  - parameter
k ∈ ℤ  - eigenvalue index or range of eigenvalues indices

"""

function charB(q::Real, k::StepRangeLen{<:Real}=1.0:1.0)
    all(x -> x > 0, k) || throw(DomainError(k, "Indices must be positive integers."))
    # Boolean indices of even and odd n values
    k_ = map(x -> round(Int, x), k);
    ie = map(iseven, k_);
    io = map(!, ie);

    b = Array{Float64}(undef, length(k));
    b[ie] = charλ(q,k)[ie];
    if q>=0
        b[io] = charλ(q,k)[io];
    else
        b[io] = charλ(abs(q), (k .+ 1))[io];
    end
    return b
end

function charB(q::Real, k::Integer)
    charB(q,k:k)
end

