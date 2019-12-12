# MathieuFunctions.jl

[![Build Status](https://travis-ci.org/BBN-Q/MathieuFunctions.jl.svg?branch=master)](https://travis-ci.org/BBN-Q/MathieuFunctions.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/35924u90n3qp7isn/branch/master?svg=true)](https://ci.appveyor.com/project/matthewware/mathieufunctions-jl/branch/master)
[![Codecov](https://codecov.io/gh/jlapeyre/MathieuFunctions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jlapeyre/MathieuFunctions.jl)
[![Coverage Status](https://coveralls.io/repos/github/BBN-Q/MathieuFunctions.jl/badge.svg?branch=master)](https://coveralls.io/github/BBN-Q/MathieuFunctions.jl?branch=master)

Julia package for Mathieu functions.

Mathieu functions are the eigenvalues and eigenfunction of *Mathieu's
equation* (for more details, see [NIST's Digital Library of
Mathematical Functions](http://dlmf.nist.gov/28)).

## Implementation

This implementation is based on 
[an algorithm described by Shirts](http://dl.acm.org/citation.cfm?id=155796), which
revolves around an application of Floquet's theorem. 

Interestingly, an
equivalent description can be obtained by considering the eigenstates
and eigenenergies of a particular type of artificial atom based on
superconducting circuits 
(see [Cotet's Ph.D. thesis](http://www.phys.ens.fr/~cottet/ACottetThesis.pdf) 
for a discussion of this connection, or 
[a recent paper by Koch et al.](https://arxiv.org/abs/cond-mat/0703002)). 
In this physical picture, Mathieu
functions can be approximated by computing the eigenvectors and
eigenvalues of a Hamiltonian truncated to a sufficiently high
dimension (this Hamiltonian corresponds to a symmetric tri-diagonal
matrix, so this approximation is rather convenient for numerics).

The current tests consist of checks against a table generated by
Mathematica (which, to our knowledge, has the most complete
implementation of Mathieu functions amongst commonly accessible
mathematical software packages).

## Mathieu equation characteristic values (eigenvalues)

We follow the notation used in NIST's Digital Library of Mathematical Functions, 

| DLMF | Julia call | Description | 
|------|------------|-------------|
| aⱼ(q)  | `charA(q,j)` | Integer order Mathieu characteristic value of even parity | 
| bⱼ(q)  | `charB(q,j)` | Integer order Mathieu characteristic value of odd parity |
| λ₍ᵥ₊ⱼ₎(q) | `charλ(q,v,j)`  | Real-valued order Mathieu characteristic value |

For integer order characteristic values, our implementation yields

![](doc/reproduction-of-dlmf-28.2.1.png?raw=true "Integer order Mathieu characteristic function")

which should be contrasted with the [corresponding plot in the DTMF @ NIST](http://dlmf.nist.gov/28.2.F1.mag).

Similarly, for real-valued order, our implementation yields

![](doc/reproduction-of-dlmf-28.13.1.png?raw=true "Real-valued order Mathieu characteristic function")

which should be contrasted with the [corresponding plot in the DTMF @ NIST](http://dlmf.nist.gov/28.13.F1.mag).

## Application: Eigenstates and energies for transmon and quantronium qubits

Thanks to the way that the Mathieu characteristic values are computed, the
order indexing convention used here is particularly convenient for computing
the energy spectrum of a broad class of superconducting qubits including the 
[quantronium](http://www.phys.ens.fr/~cottet/ACottetThesis.pdf) and the 
[transmon](https://arxiv.org/abs/cond-mat/0703002).

As an example of how this package may be used, here is the Julia code
that reproduces the energy spectrum illustrated in [Koch et al.](https://arxiv.org/abs/cond-mat/0703002)

```julia
using MathieuFunctions, PyPlot

EJoverEC = 30;
EC       = 0.35; # GHz
nLevels  = 10;
nCharge  = 2; 
nData    = 1001;

Nstart = 5;

q = -(1/2)*EJoverEC;
tol = 1e-12;

ng = LinRange(-nCharge,nCharge,nData);
nu = -2*ng;

a = zeros(nData,nLevels);
N = zeros(nData,1);
Err = zeros(nData,1);
for i=1:nData
    #println(nu[i])
    a[i,:] = charλ(q,nu[i],k=1:nLevels);
end

E = EC*a;
Eng_half = charλ(q,-1,k=1:nLevels);
Eng_0    = charλ(q, 0,k=1:nLevels);
E01      = Eng_half[2] - Eng_half[1];

plot(ng,(E .- Eng_0[1])/E01);xlabel(L"$n_g$");ylabel(L"$E$");
ylim(-.5,5)
```

![](doc/transmon-levels.png?raw=true "Transmon energy levels")

## TODO

- [ ] Some documentation with examples
- [X] Support for characteristic values (integer and non-integer order)
  - [ ] Unify order indexing, match DLMF conventions
- [ ] Support for characteristic function (integer and non-integer order)

## License

MIT (See [TLDR](https://tldrlegal.com/license/mit-license) and ["The MIT License, line by line"](https://writing.kemitchell.com/2016/09/21/MIT-License-Line-by-Line.html). )

## Contributors

Daniel Greenbaum and Marcus P da Silva

## Acknowledgements

This effort was partially supported by IARPA under the LogiQ program.

## Copyright

Raytheon BBN Technologies
