
""" 
Characteristicλ(ν,q,k)

Characteristic value λ_k for Mathieu's equation 

y'' + (λ_k - 2 q cos( 2z )) y = 0

where

ν ∈ [-1,1]  - reduced order (or characteristic exponent) of the Mathieu equation
q ∈ ℝ       - parameter
k ∈ ℤ⁺      - eigenvalue index

"""
function Characteristicλ(nu_::Real, q::Real; k=1:5) # reduced = true
    #nu = reduced ? rem(nu_+1,2)-1 : nu_;
    nu = rem(nu_+1,2)-1;
    
    # Set matrix size using formula from Shirts paper (1993), Eqs. (2.1)-(2.2).
    nu0 = nu + maximum(k);
    C = (8.46 + 0.444*nu0)/(1 + 0.085*nu0);
    D =  (0.24 + 0.0214*nu0)/(1 + 0.059*nu0);
    N = round(Int,ceil((nu0 + 2 + C*abs(q)^D)/2)); # matrix size is 2N+1
    
    d0 = (2*[-N:N;]-nu).^2;
    d1 = q*ones(2N);
    A = SymTridiagonal(d0,d1)
    a = eigvals(A,k)
    return a
end

""" 
CharacteristicA(q; k=0:4)

Characteristic value A_k for Mathieu's equation 

y'' + (A_k - 2 q cos( 2z )) y = 0

where

q ∈ ℝ  - parameter
k ∈ ℤ⁺ - eigenvalue index

"""
function CharacteristicA(q::Real; k=0:4)
    @assert all(k.>=0) "Indices must be non-negative integers."

    # Boolean indices of even and odd n values
    ie = map(iseven, k);
    io = !ie;

    a = Array(Float64,length(k));
    a[ie] = Characteristicλ(0.0,abs(q),k=k+1)[ie];
    if q>=0
        a[io] = Characteristicλ(1.0,q,k=k+1)[io]; 
    else
        if 0 in k # maybe not the cleanest way to do it
            a[io] = Characteristicλ(1.0,abs(q),k=k[2]:last(k))[io[2:end]]
        else
            a[io] = Characteristicλ(1.0,abs(q),k=k)[io]; 
        end
    end
    return a
end

""" 
CharacteristicB(q,k)

Characteristic value B_k for Mathieu's equation 

y'' + (B_k - 2 q cos( 2z )) y = 0

where

q ∈ ℝ  - parameter
k ∈ ℤ  - eigenvalueindex

"""
function CharacteristicB(q::Real; k=1:5)
    @assert all(k.>0) "Indices must be positive integers."
    # Boolean indices of even and odd n values
    ie = map(iseven, k);
    io = !ie;

    b = Array(Float64,length(k));
    b[ie] = Characteristicλ(0.0,q,k=k)[ie]; 
    if q>=0
        b[io] = Characteristicλ(1.0,q,k=k)[io];
    else
        b[io] = Characteristicλ(1.0,abs(q),k=k+1)[io]; 
    end

    return b
end


