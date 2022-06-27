---
output:
  pdf_document: default
  html_document: default
---
# Superfast response theory

$$
\newcommand{\fvec}{\mathbf{f}}
\newcommand{\evec}{\mathbf{e}}
\newcommand{\avec}{\mathbf{a}}
\newcommand{\amat}{\mathbf{A}}
\newcommand{\umat}{\mathbf{U}}
\newcommand{\kmat}{\mathbf{K}}
\newcommand{\cmat}{\mathbf{C}}
$$

## Response to edge mutations

A pair of forces acting on nodes $j$ and $k$ in opposite directions is resprsented by a vector
$$
(\fvec^{jk})^T = (\ldots, -f^{jk}\evec^{jk},\ldots,f^{jk}\evec^{jk},\ldots)
$$
Let 
$$
\avec = \amat \fvec
$$
Be a response vector ($\amat$ may be, for instance, the identity matrix, $\cmat^{1/2}$, $\cmat$, etc.). Then, the response is given by:
$$
a_i^{jk} = (\amat_{ik} - \amat_{ij})\evec_{jk}f_{jk}
$$


Analogously, the response along mode $n$ is given by: 
$$
a_n^{jk} = (\amat_{nk} - \amat_{nj})\evec_{jk}f_{jk}
$$
where $A_{ni}$ are elements of the matrix $\umat^T \amat$, where $\umat$ is the matrix of normal-mode vectors $<i|n>$ .

(Note, if $\umat$ diagonalises $\amat$ and $\alpha$ is the diagonal matrix of eigenvalues: $\umat^T\amat = \alpha \umat^T$).

If $f_{jk}$ have a distrbution with mean $0$ and standard deviation $\sigma$, then the average response is zero, and the average square response is:
$$
R_i^{jk} = <(a_i^{jk})^2> = \sigma^2 ||(\amat_{ik} - \amat_{ij})\evec_{jk}||^2
$$
and the average square response along modes:
$$
R_n^{jk} = <(a_n^{jk})^2> = \sigma^2 ||(\amat_{nk} - \amat_{nj})\evec_{jk}||^2
$$


## Response to site mutations

In the LFENM, mutating a site amounts to mutating the edges of that site. If the mutations are independent and have the same distribution of forces, then the response to a site mutation is equal to the sum of responses over edge mutations:
$$
R_i^{j} = \sum_{k \sim j} R_i^{jk} \\
R_n^{j} = \sum_{n \sim j} R_i^{jk}
$$
If the distribution of forces is not the same (for instance, in the LFENM if we use $\delta l_{jk}$ with the same distributions, $f_{jk} = k_{jk}\delta l_{jk}$ will not be identically distributed), then:
$$
R_i^{j} = \sum_{k \sim j} R_i^{jk}<f_{jk}^2> \\
R_n^{j} = \sum_{n \sim j} R_i^{jk}<f_{jk}^2>
$$


## Normal-mode transform

Note that in the normal-mode resporesentation, $\avec \rightarrow \umat^T \avec$.

## Superfast calculation

The formulas in this note allow the calculation of average response without the need to actually simulate mutations and average over them, thus in principle it should be faster (and more accurate). 

The procedure would be:

```
1. Calculate a(i, kl) and a(n, kl) for all edges (note: a(n,) = t(U) %*% a(i,))
2. Calculate r(x, kl) = a(x, kl)^2
3. calculate r(x, k) = sum(r(x, kl) * sigma_f(kl)^2)
```

