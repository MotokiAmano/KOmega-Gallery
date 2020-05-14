 ## Description

Sample code for obtaining eigenvalues using Sakurai-Sugiura (SS) method. 

## How to compile 

Using com.sh, one can compile the sample code in the supercomputer at ISSP (sekirei).
If oun wants to use other systems,
please change compile options for lapack and komega (lib and include).

By performing,
```
source ./com.sh
```
SSkomega is generated.

## Calculation
To perform calculations, 
a input file for Hamiltonian is necessary (MatrixMarket format). 
In this sample code, we used Ham.dat, which descrises
the 12-site one-dimensional Heisenberg chain.

By performing,
```
./SSKomega
```
one can perform Sakurai-Sugira method and
obtain the eigenvalues and eigenvectors
in the specified region.

## Parameters in SS method
Parameters are directly specifid in the source code
(SSKomega.c).


Parameters used Komega are following two. 
```
itermax: maximum value for iterations in komega
threshold : threshold for convergence in komega
```

The contour intergrals are performed
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

Following parameters are used in SS method.
```
nr: number of vectors used in SS method
nk: order  of moment used in  SS method
```

Author: Takahiro Misawa (ISSP, Univ. of Tokyo), Date: 2020/1/8
