 ## Description

Sample code for perfoming BiCG method. 

## How to compile 

Using com.sh, one can compile the sample code in the supercomputer at ISSP (sekirei).
If oun wants to use other systems,
please change compile options for lapack and komega (lib and include).

By performing,
```
source ./com.sh
```
Komega.out is generated.

## Calculation
To perform calculations, 
a input file for Hamiltonian is necessary (MatrixMarket format). 
In this sample code, we used Ham.dat, which descrises
the 12-site one-dimensional Heisenberg chain.

By performing,
```
./Komega.out
```
one can perform BiCG method and
obtain the solutions of (zI-H)x=b.

## Parameters in BiCG method
Parameters are directly specifid in the source code
(Komega.c).


Parameters used Komega are following two. 
```
itermax: maximum value for iterations in komega
threshold : threshold for convergence in komega
```

Complex numbers z are selected
at the following points.

```
zj=γ+ρ*exp[(2πI/nz)*(j+1/2)]
```

Here, 
```
gamma: origin of the contour  integrals
rho: length of the contour  integrals
nz:  points of the contour intergral
```

Author: Takahiro Misawa (ISSP, Univ. of Tokyo), Date: 2020/5/14
